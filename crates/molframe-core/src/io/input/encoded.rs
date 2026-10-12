//! Bounded decoding from an existing allocation: O(expanded bytes), no encoded copy.

use super::InputBuffer;
use super::collect::OutputLimit;
#[cfg(not(feature = "mmap"))]
use super::decode::decompress;
#[cfg(all(feature = "mmap", feature = "gzip"))]
use super::decode::map_gzip_reader;
#[cfg(all(feature = "mmap", feature = "zstd"))]
use super::decode::map_zstd_reader;
use super::{Compression, Limits};
use crate::diagnostic::Diagnostic;

impl InputBuffer {
    /// Adopts plain bytes or decodes a gzip/Zstandard container within `limits`.
    ///
    /// Plain bytes retain their allocation. Compressed bytes are decoded directly
    /// from the supplied vector; expanded output uses private mapped storage when
    /// `mmap` is enabled. [`Self::from_bytes`] deliberately performs no decoding.
    ///
    /// ```
    /// use molframe_core::io::{InputBuffer, Limits};
    /// let bytes = b"data_example\n".to_vec();
    /// let pointer = bytes.as_ptr();
    /// let input = InputBuffer::from_encoded_bytes(bytes, Limits::default())?;
    /// assert_eq!(input.as_bytes().as_ptr(), pointer);
    /// # Ok::<(), molframe_core::diagnostic::Diagnostic>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Invalid containers, unavailable decoders, or exceeded byte/expansion limits.
    pub fn from_encoded_bytes(bytes: Vec<u8>, limits: Limits) -> Result<Self, Diagnostic> {
        let length = u64::try_from(bytes.len())
            .map_err(|_| Limits::exceeded("decompressed bytes", bytes.len()))?;
        OutputLimit::uncompressed(limits).validate(length)?;
        let compression = Compression::sniff(&bytes);
        if compression == Compression::None {
            return Ok(Self::from_bytes(bytes));
        }
        #[cfg(not(feature = "mmap"))]
        {
            Ok(Self::from_bytes(decompress(bytes, limits)?))
        }
        #[cfg(feature = "mmap")]
        {
            match compression {
                #[cfg(feature = "gzip")]
                Compression::Gzip => {
                    map_gzip_reader(bytes.as_slice(), length, limits, None).map(Self::from_mapped)
                }
                #[cfg(not(feature = "gzip"))]
                Compression::Gzip => Err(super::decode::unsupported("gzip")),
                #[cfg(feature = "zstd")]
                Compression::Zstd => {
                    map_zstd_reader(bytes.as_slice(), length, limits, None).map(Self::from_mapped)
                }
                #[cfg(not(feature = "zstd"))]
                Compression::Zstd => Err(super::decode::unsupported("zstd")),
                Compression::None => Ok(Self::from_bytes(bytes)),
            }
        }
    }
}

#[cfg(test)]
#[path = "encoded_tests.rs"]
mod tests;
