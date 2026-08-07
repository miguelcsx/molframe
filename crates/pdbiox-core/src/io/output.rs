//! Deterministic output encoding selected by the destination suffix.

use crate::diagnostic::{Code, Diagnostic};
use std::path::Path;

/// Writes bytes to a destination, applying a named compression container.
///
/// Gzip uses a zero timestamp, so identical input produces identical bytes.
///
/// # Errors
///
/// Returns `E7901` when compression is unavailable or the destination cannot
/// be written.
pub fn write_output(path: impl AsRef<Path>, bytes: &[u8]) -> Result<(), Diagnostic> {
    let path = path.as_ref();
    let encoded = encode(path, bytes)?;
    std::fs::write(path, encoded).map_err(|error| output_error(path, error))
}

fn encode(path: &Path, bytes: &[u8]) -> Result<Vec<u8>, Diagnostic> {
    match path.extension().and_then(|suffix| suffix.to_str()) {
        Some(suffix) if suffix.eq_ignore_ascii_case("gz") => gzip(path, bytes),
        Some(suffix) if suffix.eq_ignore_ascii_case("zst") => zstd(path, bytes),
        _ => Ok(bytes.to_vec()),
    }
}

#[cfg(feature = "gzip")]
fn gzip(path: &Path, bytes: &[u8]) -> Result<Vec<u8>, Diagnostic> {
    use flate2::{Compression, GzBuilder};
    use std::io::Write;

    let mut encoder = GzBuilder::new()
        .mtime(0)
        .write(Vec::new(), Compression::default());
    encoder
        .write_all(bytes)
        .map_err(|error| output_error(path, error))?;
    encoder.finish().map_err(|error| output_error(path, error))
}

#[cfg(not(feature = "gzip"))]
fn gzip(path: &Path, _bytes: &[u8]) -> Result<Vec<u8>, Diagnostic> {
    Err(unavailable(path, "gzip"))
}

#[cfg(feature = "zstd")]
fn zstd(path: &Path, bytes: &[u8]) -> Result<Vec<u8>, Diagnostic> {
    zstd::stream::encode_all(bytes, 0).map_err(|error| output_error(path, error))
}

#[cfg(not(feature = "zstd"))]
fn zstd(path: &Path, _bytes: &[u8]) -> Result<Vec<u8>, Diagnostic> {
    Err(unavailable(path, "zstandard"))
}

#[cfg(any(not(feature = "gzip"), not(feature = "zstd")))]
fn unavailable(path: &Path, compression: &str) -> Diagnostic {
    Diagnostic::new(Code::E7901)
        .with_context("path", path.display().to_string())
        .with_context("compression", compression)
        .with_context("reason", "support is not enabled in this build")
}

fn output_error(path: &Path, error: impl std::fmt::Display) -> Diagnostic {
    Diagnostic::new(Code::E7901)
        .with_context("path", path.display().to_string())
        .with_context("reason", error.to_string())
}

#[cfg(test)]
#[path = "output_tests.rs"]
mod tests;
