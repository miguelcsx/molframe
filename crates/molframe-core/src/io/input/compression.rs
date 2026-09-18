//! Compression-container recognition by magic bytes.

/// How a stream of bytes was compressed.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[non_exhaustive]
pub enum Compression {
    /// Not compressed.
    None,
    /// Gzip.
    Gzip,
    /// Zstandard.
    Zstd,
}

impl Compression {
    /// Recognises a container without consulting its file name.
    #[must_use]
    #[inline]
    pub fn sniff(bytes: &[u8]) -> Self {
        match bytes {
            [0x1f, 0x8b, ..] => Self::Gzip,
            [0x28, 0xb5, 0x2f, 0xfd, ..] => Self::Zstd,
            _ => Self::None,
        }
    }
}
