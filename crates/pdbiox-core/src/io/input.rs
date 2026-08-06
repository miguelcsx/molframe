//! Where bytes come from, and the limits that keep a hostile file from
//! exhausting the machine.
//!
//! Memory mapping is right for a large uncompressed file on local disk and wrong
//! for a compressed one, a network stream or bytes already in hand. All four are
//! ordinary here and none is architecturally privileged: a reader takes bytes and
//! does not care how they were obtained.
//!
//! Decompression is transparent and bounded. A file that expands a thousandfold
//! is refused rather than obeyed.

use crate::diagnostic::{Code, Diagnostic};
use std::fs::File;
use std::io::Read;
use std::path::Path;
use std::sync::Arc;

/// Ceilings that apply while reading.
///
/// These exist to turn a hostile input into a diagnostic rather than an
/// out-of-memory kill. The defaults are far above any real file.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Limits {
    /// Largest accepted size after decompression, in bytes.
    pub decompressed_bytes: u64,
    /// Largest accepted ratio of decompressed to compressed size.
    pub compression_ratio: u64,
    /// Largest accepted declared row count for one category.
    pub rows_per_category: u64,
    /// Deepest accepted nesting.
    pub nesting_depth: u32,
    /// Largest accepted identifier dictionary.
    pub dictionary_entries: u32,
}

impl Default for Limits {
    fn default() -> Self {
        Self {
            decompressed_bytes: 4 << 30,
            compression_ratio: 1_000,
            rows_per_category: 100_000_000,
            nesting_depth: 64,
            dictionary_entries: 1_000_000,
        }
    }
}

impl Limits {
    /// The finding raised when a limit is exceeded.
    #[must_use]
    pub fn exceeded(limit: &'static str, value: impl std::fmt::Display) -> Diagnostic {
        Diagnostic::new(Code::E1901)
            .with_context("limit", limit)
            .with_context("value", value.to_string())
    }
}

/// How a stream of bytes was compressed, where it was.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[non_exhaustive]
pub enum Compression {
    /// Not compressed.
    None,
    /// Gzip, as in a `.gz` suffix.
    Gzip,
    /// Zstandard, as in a `.zst` suffix.
    Zstd,
}

impl Compression {
    /// Recognises a compression container from its leading bytes.
    ///
    /// Content decides, not the file name: an archive mirror serving gzip under
    /// a plain `.cif` name is common enough to matter.
    #[must_use]
    pub fn sniff(bytes: &[u8]) -> Self {
        match bytes {
            [0x1f, 0x8b, ..] => Self::Gzip,
            [0x28, 0xb5, 0x2f, 0xfd, ..] => Self::Zstd,
            _ => Self::None,
        }
    }
}

/// Bytes to read, however they were obtained.
///
/// # Examples
///
/// ```
/// use pdbiox_core::io::InputBuffer;
///
/// let input = InputBuffer::from_bytes(b"data_test\n".to_vec());
/// assert_eq!(input.len(), 10);
/// assert!(input.as_bytes().starts_with(b"data_"));
/// ```
#[derive(Clone, Debug)]
pub struct InputBuffer {
    bytes: Arc<[u8]>,
    origin: Option<Arc<str>>,
}

impl InputBuffer {
    /// Takes ownership of bytes already in hand.
    #[must_use]
    pub fn from_bytes(bytes: Vec<u8>) -> Self {
        Self {
            bytes: bytes.into(),
            origin: None,
        }
    }

    /// Reads a file, decompressing it if it is compressed.
    ///
    /// # Errors
    ///
    /// Returns a finding when the file cannot be read, or when decompressing it
    /// would exceed `limits`.
    pub fn open(path: impl AsRef<Path>, limits: Limits) -> Result<Self, Diagnostic> {
        let path = path.as_ref();
        let mut file = File::open(path).map_err(|error| {
            Diagnostic::new(Code::E1901)
                .with_message("input could not be opened")
                .with_context("path", path.display().to_string())
                .with_context("reason", error.to_string())
        })?;

        let mut raw = Vec::new();
        file.read_to_end(&mut raw).map_err(|error| {
            Diagnostic::new(Code::E1901)
                .with_message("input could not be read")
                .with_context("path", path.display().to_string())
                .with_context("reason", error.to_string())
        })?;

        let bytes = decompress(raw, limits)?;
        Ok(Self {
            bytes: bytes.into(),
            origin: Some(Arc::from(path.display().to_string().as_str())),
        })
    }

    /// The bytes.
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        &self.bytes
    }

    /// The number of bytes.
    #[must_use]
    pub fn len(&self) -> usize {
        self.bytes.len()
    }

    /// Returns true when there is nothing to read.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.bytes.is_empty()
    }

    /// Where the bytes came from, for a diagnostic to name.
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

/// Expands a compressed buffer, refusing one that expands beyond the limits.
fn decompress(raw: Vec<u8>, limits: Limits) -> Result<Vec<u8>, Diagnostic> {
    let compressed = raw.len() as u64;
    match Compression::sniff(&raw) {
        Compression::None => Ok(raw),
        Compression::Gzip => expand_gzip(&raw, compressed, limits),
        Compression::Zstd => expand_zstd(&raw, compressed, limits),
    }
}

#[cfg(feature = "gzip")]
fn expand_gzip(raw: &[u8], compressed: u64, limits: Limits) -> Result<Vec<u8>, Diagnostic> {
    use std::io::Read as _;
    let mut decoder = flate2::read::GzDecoder::new(raw);
    let mut out = Vec::new();
    decoder
        .by_ref()
        .take(limits.decompressed_bytes)
        .read_to_end(&mut out)
        .map_err(|error| {
            Diagnostic::new(Code::E1901)
                .with_message("gzip stream could not be decoded")
                .with_context("reason", error.to_string())
        })?;
    check_expansion(out.len() as u64, compressed, limits)?;
    Ok(out)
}

#[cfg(not(feature = "gzip"))]
fn expand_gzip(_: &[u8], _: u64, _: Limits) -> Result<Vec<u8>, Diagnostic> {
    Err(unsupported("gzip"))
}

#[cfg(feature = "zstd")]
fn expand_zstd(raw: &[u8], compressed: u64, limits: Limits) -> Result<Vec<u8>, Diagnostic> {
    let out = zstd::decode_all(raw).map_err(|error| {
        Diagnostic::new(Code::E1901)
            .with_message("zstd stream could not be decoded")
            .with_context("reason", error.to_string())
    })?;
    check_expansion(out.len() as u64, compressed, limits)?;
    Ok(out)
}

#[cfg(not(feature = "zstd"))]
fn expand_zstd(_: &[u8], _: u64, _: Limits) -> Result<Vec<u8>, Diagnostic> {
    Err(unsupported("zstd"))
}

#[cfg(any(not(feature = "gzip"), not(feature = "zstd")))]
fn unsupported(container: &'static str) -> Diagnostic {
    Diagnostic::new(Code::E1901)
        .with_message("input is compressed with a container this build does not include")
        .with_context("container", container)
}

#[cfg_attr(
    not(any(feature = "gzip", feature = "zstd")),
    expect(dead_code, reason = "only a compression feature reaches this")
)]
fn check_expansion(expanded: u64, compressed: u64, limits: Limits) -> Result<(), Diagnostic> {
    if expanded > limits.decompressed_bytes {
        return Err(Limits::exceeded("decompressed bytes", expanded));
    }
    if compressed > 0 && expanded / compressed.max(1) > limits.compression_ratio {
        return Err(Limits::exceeded(
            "compression ratio",
            expanded / compressed.max(1),
        ));
    }
    Ok(())
}

#[cfg(test)]
#[path = "input_tests.rs"]
mod tests;
