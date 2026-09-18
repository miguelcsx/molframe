//! Compression-prefix replay and bounded gzip/Zstandard decoding.

#[cfg(all(feature = "mmap", any(feature = "gzip", feature = "zstd")))]
use super::collect::map_checked;
use super::collect::read_failure;
#[cfg(any(feature = "gzip", feature = "zstd"))]
use super::collect::{OutputLimit, collect_checked};
use super::{Compression, Limits};
use crate::diagnostic::{Code, Diagnostic};
use std::io::{ErrorKind, Read};
use std::path::Path;

const COMPRESSION_PREFIX_BYTES: usize = 4;

pub(super) struct ReplayReader<R> {
    prefix: [u8; COMPRESSION_PREFIX_BYTES],
    prefix_len: usize,
    position: usize,
    inner: R,
}

impl<R> ReplayReader<R> {
    const fn new(prefix: [u8; COMPRESSION_PREFIX_BYTES], prefix_len: usize, inner: R) -> Self {
        Self {
            prefix,
            prefix_len,
            position: 0,
            inner,
        }
    }
}

impl<R: Read> Read for ReplayReader<R> {
    fn read(&mut self, buffer: &mut [u8]) -> std::io::Result<usize> {
        if self.position < self.prefix_len {
            let count = (self.prefix_len - self.position).min(buffer.len());
            if count == 0 {
                return Ok(0);
            }
            let end = self.position + count;
            buffer[..count].copy_from_slice(&self.prefix[self.position..end]);
            self.position = end;
            return Ok(count);
        }
        self.inner.read(buffer)
    }
}

pub(super) fn sniff_reader<R: Read>(
    mut reader: R,
    path: &Path,
) -> Result<(Compression, ReplayReader<R>), Diagnostic> {
    let mut prefix = [0_u8; COMPRESSION_PREFIX_BYTES];
    let mut prefix_len = 0;
    while prefix_len < prefix.len() {
        match reader.read(&mut prefix[prefix_len..]) {
            Ok(0) => break,
            Ok(count) => prefix_len += count,
            Err(error) if error.kind() == ErrorKind::Interrupted => {}
            Err(error) => return Err(read_failure("input could not be read", Some(path), &error)),
        }
    }
    Ok((
        Compression::sniff(&prefix[..prefix_len]),
        ReplayReader::new(prefix, prefix_len, reader),
    ))
}

pub(super) fn compressed_file_size(
    size: Result<u64, std::io::Error>,
    path: &Path,
) -> Result<u64, Diagnostic> {
    let size =
        size.map_err(|error| read_failure("input size could not be read", Some(path), &error))?;
    if size == 0 {
        return Err(Diagnostic::new(Code::E1901)
            .with_message("compressed input size could not be determined")
            .with_context("path", path.display().to_string()));
    }
    Ok(size)
}

#[cfg(not(feature = "mmap"))]
pub(super) fn read_bounded(mut reader: impl Read, limits: Limits) -> std::io::Result<Vec<u8>> {
    let mut raw = Vec::new();
    match limits.decompressed_bytes.checked_add(1) {
        Some(ceiling) => {
            reader.take(ceiling).read_to_end(&mut raw)?;
        }
        None => {
            reader.read_to_end(&mut raw)?;
        }
    }
    let length = u64::try_from(raw.len())
        .map_err(|_| std::io::Error::other("input length exceeds u64::MAX"))?;
    if length > limits.decompressed_bytes {
        return Err(std::io::Error::other(
            "input exceeds the configured byte limit",
        ));
    }
    Ok(raw)
}

#[cfg(not(feature = "mmap"))]
pub(super) fn decompress(raw: Vec<u8>, limits: Limits) -> Result<Vec<u8>, Diagnostic> {
    let Ok(compressed) = u64::try_from(raw.len()) else {
        return Err(Limits::exceeded("decompressed bytes", raw.len()));
    };
    match Compression::sniff(&raw) {
        Compression::None => Ok(raw),
        Compression::Gzip => expand_gzip(&raw, compressed, limits),
        Compression::Zstd => expand_zstd(&raw, compressed, limits),
    }
}

#[cfg(all(feature = "gzip", not(feature = "mmap")))]
fn expand_gzip(raw: &[u8], compressed: u64, limits: Limits) -> Result<Vec<u8>, Diagnostic> {
    expand_gzip_reader(raw, compressed, limits, None)
}

#[cfg(all(not(feature = "gzip"), not(feature = "mmap")))]
fn expand_gzip(_: &[u8], _: u64, _: Limits) -> Result<Vec<u8>, Diagnostic> {
    Err(unsupported("gzip"))
}

#[cfg(feature = "gzip")]
pub(super) fn expand_gzip_reader<R: Read>(
    reader: R,
    compressed: u64,
    limits: Limits,
    path: Option<&Path>,
) -> Result<Vec<u8>, Diagnostic> {
    collect_checked(
        flate2::read::GzDecoder::new(reader),
        OutputLimit::expanded(compressed, limits),
        "gzip stream could not be decoded",
        path,
    )
}

#[cfg(all(feature = "gzip", feature = "mmap"))]
pub(super) fn map_gzip_reader<R: Read>(
    reader: R,
    compressed: u64,
    limits: Limits,
    path: Option<&Path>,
) -> Result<molframe_mmap::MappedFile, Diagnostic> {
    map_checked(
        flate2::read::GzDecoder::new(reader),
        OutputLimit::expanded(compressed, limits),
        "gzip stream could not be decoded",
        path,
    )
}

#[cfg(not(feature = "gzip"))]
pub(super) fn expand_gzip_reader<R: Read>(
    _: R,
    _: u64,
    _: Limits,
    _: Option<&Path>,
) -> Result<Vec<u8>, Diagnostic> {
    Err(unsupported("gzip"))
}

#[cfg(all(feature = "zstd", not(feature = "mmap")))]
fn expand_zstd(raw: &[u8], compressed: u64, limits: Limits) -> Result<Vec<u8>, Diagnostic> {
    expand_zstd_reader(raw, compressed, limits, None)
}

#[cfg(all(not(feature = "zstd"), not(feature = "mmap")))]
fn expand_zstd(_: &[u8], _: u64, _: Limits) -> Result<Vec<u8>, Diagnostic> {
    Err(unsupported("zstd"))
}

#[cfg(feature = "zstd")]
pub(super) fn expand_zstd_reader<R: Read>(
    reader: R,
    compressed: u64,
    limits: Limits,
    path: Option<&Path>,
) -> Result<Vec<u8>, Diagnostic> {
    let decoder = zstd::stream::read::Decoder::new(reader)
        .map_err(|error| read_failure("zstd stream could not be decoded", path, &error))?;
    collect_checked(
        decoder,
        OutputLimit::expanded(compressed, limits),
        "zstd stream could not be decoded",
        path,
    )
}

#[cfg(all(feature = "zstd", feature = "mmap"))]
pub(super) fn map_zstd_reader<R: Read>(
    reader: R,
    compressed: u64,
    limits: Limits,
    path: Option<&Path>,
) -> Result<molframe_mmap::MappedFile, Diagnostic> {
    let decoder = zstd::stream::read::Decoder::new(reader)
        .map_err(|error| read_failure("zstd stream could not be decoded", path, &error))?;
    map_checked(
        decoder,
        OutputLimit::expanded(compressed, limits),
        "zstd stream could not be decoded",
        path,
    )
}

#[cfg(not(feature = "zstd"))]
pub(super) fn expand_zstd_reader<R: Read>(
    _: R,
    _: u64,
    _: Limits,
    _: Option<&Path>,
) -> Result<Vec<u8>, Diagnostic> {
    Err(unsupported("zstd"))
}

#[cfg(any(not(feature = "gzip"), not(feature = "zstd")))]
pub(super) fn unsupported(container: &'static str) -> Diagnostic {
    Diagnostic::new(Code::E1901)
        .with_message("input is compressed with a container this build does not include")
        .with_context("container", container)
}
