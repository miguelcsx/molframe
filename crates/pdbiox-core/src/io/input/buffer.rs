//! Shared storage and public input construction.

use super::collect::{OutputLimit, collect_checked, path_origin, read_failure};
use super::decode::{
    compressed_file_size, decompress, expand_gzip_reader, expand_zstd_reader, read_bounded,
    sniff_reader,
};
use super::{Compression, Limits};
use crate::diagnostic::{Code, Diagnostic};
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
    /// Read-only pages mapped from a local file.
    #[cfg(feature = "mmap")]
    Mapped,
}

#[derive(Clone, Debug)]
enum InputBytes {
    Owned(Arc<[u8]>),
    #[cfg(feature = "mmap")]
    Mapped(Arc<pdbiox_mmap::MappedFile>),
}

impl InputBytes {
    fn as_slice(&self) -> &[u8] {
        match self {
            Self::Owned(bytes) => bytes,
            #[cfg(feature = "mmap")]
            Self::Mapped(mapped) => mapped.as_bytes(),
        }
    }

    const fn kind(&self) -> InputKind {
        match self {
            Self::Owned(_) => InputKind::Owned,
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
            bytes: InputBytes::Owned(bytes.into()),
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
        if let Ok(length) = compressed_size.as_ref()
            && *length <= limits.decompressed_bytes
            && *length >= MMAP_MIN_BYTES
            && let Ok(mapped) = pdbiox_mmap::MappedFile::new(&file)
            && Compression::sniff(mapped.as_bytes()) == Compression::None
        {
            return Ok(Self {
                bytes: InputBytes::Mapped(Arc::new(mapped)),
                origin: Some(path_origin(path)),
            });
        }

        let (compression, reader) = sniff_reader(file, path)?;
        let compressed = match compression {
            Compression::None => 0,
            Compression::Gzip | Compression::Zstd => compressed_file_size(compressed_size, path)?,
        };
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
            bytes: InputBytes::Owned(bytes.into()),
            origin: Some(path_origin(path)),
        })
    }

    /// Reads a bounded stream and transparently decompresses its bytes.
    ///
    /// # Errors
    ///
    /// Returns a registered resource or limit diagnostic.
    pub fn from_reader(reader: impl Read, limits: Limits) -> Result<Self, Diagnostic> {
        let raw = read_bounded(reader, limits).map_err(|error| {
            Diagnostic::new(Code::E1901)
                .with_message("input stream could not be read")
                .with_context("reason", error.to_string())
        })?;
        Ok(Self {
            bytes: InputBytes::Owned(decompress(raw, limits)?.into()),
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
