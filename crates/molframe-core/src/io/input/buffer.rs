//! Shared storage and public input construction.

#[cfg(feature = "mmap")]
use super::collect::map_checked;
use super::collect::{OutputLimit, collect_checked, path_origin, read_failure};
#[cfg(all(feature = "gzip", feature = "mmap"))]
use super::decode::map_gzip_reader;
#[cfg(all(feature = "zstd", feature = "mmap"))]
use super::decode::map_zstd_reader;
use super::decode::{compressed_file_size, expand_gzip_reader, expand_zstd_reader, sniff_reader};
#[cfg(not(feature = "mmap"))]
use super::decode::{decompress, read_bounded};
use super::{Compression, Limits};
#[cfg(not(feature = "mmap"))]
use crate::diagnostic::Code;
use crate::diagnostic::Diagnostic;
use std::fmt;
use std::fs::File;
use std::io::Read;
use std::path::Path;
use std::sync::Arc;

#[cfg(feature = "mmap")]
pub(super) const MMAP_MIN_BYTES: u64 = 1 << 20;

/// The storage that backs an input byte slice.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum InputKind {
    /// Bytes owned by the process.
    Owned,
    /// Stable bytes retained through an external owner without copying.
    OwnerBacked,
    /// Read-only pages mapped from a local file.
    #[cfg(feature = "mmap")]
    Mapped,
}

trait ByteOwner: AsRef<[u8]> + fmt::Debug + Send + Sync {}

impl<T> ByteOwner for T where T: AsRef<[u8]> + fmt::Debug + Send + Sync {}

#[derive(Clone, Debug)]
enum InputBytes {
    Owned(Arc<Vec<u8>>),
    OwnerBacked(Arc<dyn ByteOwner>),
    #[cfg(feature = "mmap")]
    Mapped(Arc<molframe_mmap::MappedFile>),
}

impl InputBytes {
    fn as_slice(&self) -> &[u8] {
        match self {
            Self::Owned(bytes) => bytes.as_slice(),
            Self::OwnerBacked(owner) => owner.as_ref().as_ref(),
            #[cfg(feature = "mmap")]
            Self::Mapped(mapped) => mapped.as_bytes(),
        }
    }

    const fn kind(&self) -> InputKind {
        match self {
            Self::Owned(_) => InputKind::Owned,
            Self::OwnerBacked(_) => InputKind::OwnerBacked,
            #[cfg(feature = "mmap")]
            Self::Mapped(_) => InputKind::Mapped,
        }
    }
}

/// Bytes to read, however they were obtained.
#[derive(Clone, Debug)]
pub struct InputBuffer {
    bytes: InputBytes,
    origin: Option<Arc<str>>,
}

impl InputBuffer {
    /// Takes ownership of bytes already in hand.
    #[must_use]
    pub fn from_bytes(bytes: Vec<u8>) -> Self {
        Self {
            bytes: InputBytes::Owned(Arc::new(bytes)),
            origin: None,
        }
    }

    /// Retains stable externally owned bytes without copying their payload.
    ///
    /// The owner is held for the complete lifetime of this buffer and every
    /// clone. This is the adoption path for language-runtime byte objects and
    /// existing immutable buffers whose address remains stable.
    #[must_use]
    pub fn from_owner<T>(owner: T) -> Self
    where
        T: AsRef<[u8]> + fmt::Debug + Send + Sync + 'static,
    {
        Self {
            bytes: InputBytes::OwnerBacked(Arc::new(owner)),
            origin: None,
        }
    }

    /// Opens and transparently decompresses a bounded local file.
    ///
    /// # Errors
    ///
    /// Returns a registered resource or limit diagnostic.
    pub fn open(path: impl AsRef<Path>, limits: Limits) -> Result<Self, Diagnostic> {
        let path = path.as_ref();
        let file = File::open(path)
            .map_err(|error| read_failure("input could not be opened", Some(path), &error))?;
        let compressed_size = file.metadata().map(|metadata| metadata.len());
        #[cfg(feature = "mmap")]
        let file_length = compressed_size.as_ref().copied().ok();

        let (compression, reader) = sniff_reader(file, path)?;
        let compressed = match compression {
            Compression::None => 0,
            Compression::Gzip | Compression::Zstd => compressed_file_size(compressed_size, path)?,
        };
        #[cfg(feature = "mmap")]
        if compression != Compression::None
            || file_length.is_some_and(|length| length >= MMAP_MIN_BYTES)
        {
            let mapped = match compression {
                Compression::None => map_checked(
                    reader,
                    OutputLimit::uncompressed(limits),
                    "input could not be read",
                    Some(path),
                )?,
                #[cfg(feature = "gzip")]
                Compression::Gzip => map_gzip_reader(reader, compressed, limits, Some(path))?,
                #[cfg(not(feature = "gzip"))]
                Compression::Gzip => return Err(super::decode::unsupported("gzip")),
                #[cfg(feature = "zstd")]
                Compression::Zstd => map_zstd_reader(reader, compressed, limits, Some(path))?,
                #[cfg(not(feature = "zstd"))]
                Compression::Zstd => return Err(super::decode::unsupported("zstd")),
            };
            return Ok(Self {
                bytes: InputBytes::Mapped(Arc::new(mapped)),
                origin: Some(path_origin(path)),
            });
        }

        let bytes = match compression {
            Compression::None => collect_checked(
                reader,
                OutputLimit::uncompressed(limits),
                "input could not be read",
                Some(path),
            )?,
            Compression::Gzip => expand_gzip_reader(reader, compressed, limits, Some(path))?,
            Compression::Zstd => expand_zstd_reader(reader, compressed, limits, Some(path))?,
        };
        Ok(Self {
            bytes: InputBytes::Owned(Arc::new(bytes)),
            origin: Some(path_origin(path)),
        })
    }

    /// Wraps a mapping the caller has already established.
    ///
    /// [`Self::open`] streams a source through a bounded buffer into private
    /// storage before mapping it, so opening a hundred-gigabyte entry reads a
    /// hundred gigabytes, writes a hundred gigabytes of scratch, and only then
    /// begins. A caller who can guarantee the file will not change underneath
    /// the mapping can avoid all of that by mapping it in place and handing the
    /// result here.
    ///
    /// Establishing such a mapping requires `unsafe`, because no operating
    /// system can enforce the guarantee. This crate forbids `unsafe` and so
    /// cannot offer it: the call belongs to
    /// [`molframe_mmap::MappedFile::map_file_unchecked`], in the crate that
    /// audits it, written by the caller who is in a position to make the
    /// promise. By the time a [`MappedFile`](molframe_mmap::MappedFile) exists
    /// the obligation is discharged, which is why this constructor is safe.
    ///
    /// The mapping is advised as a single forward pass, which is how every
    /// parser above this crate reads. A platform that declines the hint is not
    /// an error.
    #[cfg(feature = "mmap")]
    #[must_use]
    pub fn from_mapped(mapped: molframe_mmap::MappedFile) -> Self {
        let _advised = mapped.advise_sequential();
        Self {
            bytes: InputBytes::Mapped(Arc::new(mapped)),
            origin: None,
        }
    }

    /// Reads a bounded stream and transparently decompresses its bytes.
    ///
    /// # Errors
    ///
    /// Returns a registered resource or limit diagnostic.
    pub fn from_reader(reader: impl Read, limits: Limits) -> Result<Self, Diagnostic> {
        #[cfg(feature = "mmap")]
        {
            Self::map_reader(reader, limits)
        }
        #[cfg(not(feature = "mmap"))]
        {
            let raw = read_bounded(reader, limits).map_err(|error| {
                Diagnostic::new(Code::E1901)
                    .with_message("input stream could not be read")
                    .with_context("reason", error.to_string())
            })?;
            let bytes = decompress(raw, limits)?;
            Ok(Self {
                bytes: InputBytes::Owned(Arc::new(bytes)),
                origin: None,
            })
        }
    }

    #[cfg(feature = "mmap")]
    fn map_reader(reader: impl Read, limits: Limits) -> Result<Self, Diagnostic> {
        let raw = map_checked(
            reader,
            OutputLimit::uncompressed(limits),
            "input stream could not be read",
            None,
        )?;
        #[cfg(any(feature = "gzip", feature = "zstd"))]
        let compressed = u64::try_from(raw.as_bytes().len())
            .map_err(|_| Limits::exceeded("compressed bytes", "more than u64::MAX"))?;
        let bytes = match Compression::sniff(raw.as_bytes()) {
            Compression::None => raw,
            #[cfg(feature = "gzip")]
            Compression::Gzip => map_gzip_reader(raw.as_bytes(), compressed, limits, None)?,
            #[cfg(not(feature = "gzip"))]
            Compression::Gzip => return Err(super::decode::unsupported("gzip")),
            #[cfg(feature = "zstd")]
            Compression::Zstd => map_zstd_reader(raw.as_bytes(), compressed, limits, None)?,
            #[cfg(not(feature = "zstd"))]
            Compression::Zstd => return Err(super::decode::unsupported("zstd")),
        };
        Ok(Self {
            bytes: InputBytes::Mapped(Arc::new(bytes)),
            origin: None,
        })
    }

    /// The bytes.
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        self.bytes.as_slice()
    }

    /// The number of bytes.
    #[must_use]
    pub fn len(&self) -> usize {
        self.bytes.as_slice().len()
    }

    /// Returns true when there is nothing to read.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.bytes.as_slice().is_empty()
    }

    /// How this buffer owns its bytes.
    #[must_use]
    pub const fn kind(&self) -> InputKind {
        self.bytes.kind()
    }

    /// Where the bytes came from.
    #[must_use]
    pub fn origin(&self) -> Option<&str> {
        self.origin.as_deref()
    }

    /// Names where the bytes came from.
    #[must_use]
    pub fn with_origin(mut self, origin: &str) -> Self {
        self.origin = Some(Arc::from(origin));
        self
    }
}
