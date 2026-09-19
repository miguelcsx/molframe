//! Plain or explicitly spooled compressed files behind one window contract.

use super::{ByteWindow, Compression, SourceBytes, SpillWindowedFile, WindowedFile};
use crate::{Code, Diagnostic, ExecutionContext};
use std::fs::File;
#[cfg(any(feature = "gzip", feature = "zstd"))]
use std::io::ErrorKind;
use std::io::Read;
use std::path::Path;
#[cfg(any(feature = "gzip", feature = "zstd"))]
use std::sync::atomic::{AtomicU64, Ordering};

#[cfg(any(feature = "gzip", feature = "zstd"))]
static SPOOL_SEQUENCE: AtomicU64 = AtomicU64::new(0);

/// A seekable bounded source, with compressed input spooled only when the
/// execution context explicitly permits temporary storage.
#[derive(Debug)]
pub enum WindowedSourceFile {
    /// Uncompressed source read directly from its file.
    Plain(WindowedFile),
    /// Decompressed bytes held in execution-owned spill storage.
    Spilled(SpillWindowedFile),
}

impl WindowedSourceFile {
    /// Opens a plain, gzip or Zstandard file with bounded resident memory.
    ///
    /// Compressed inputs require an enabled [`crate::TempStoragePolicy`]. The
    /// spool is charged, checksummed and deleted with this source.
    ///
    /// # Errors
    ///
    /// Returns a registered I/O, codec, spill-budget or memory diagnostic.
    pub fn open(
        path: impl AsRef<Path>,
        max_window_bytes: usize,
        context: &ExecutionContext,
    ) -> Result<Self, Diagnostic> {
        let path = path.as_ref();
        let mut file = File::open(path).map_err(|error| io_error("open", &error))?;
        let mut prefix = [0_u8; 4];
        let prefix_len = file
            .read(&mut prefix)
            .map_err(|error| io_error("read compression prefix", &error))?;
        match Compression::sniff(&prefix[..prefix_len]) {
            Compression::None => {
                WindowedFile::open(path, max_window_bytes, context).map(Self::Plain)
            }
            Compression::Gzip => open_gzip(path, max_window_bytes, context).map(Self::Spilled),
            Compression::Zstd => open_zstd(path, max_window_bytes, context).map(Self::Spilled),
        }
    }

    /// Physical bytes currently retained in spill storage.
    #[must_use]
    pub const fn spill_bytes(&self) -> u64 {
        match self {
            Self::Plain(_) => 0,
            Self::Spilled(source) => source.spill_bytes(),
        }
    }
}

impl SourceBytes for WindowedSourceFile {
    fn window(&mut self, start: u64, len: usize) -> Result<ByteWindow<'_>, Diagnostic> {
        match self {
            Self::Plain(source) => source.window(start, len),
            Self::Spilled(source) => source.window(start, len),
        }
    }

    fn len_hint(&self) -> Option<u64> {
        match self {
            Self::Plain(source) => source.len_hint(),
            Self::Spilled(source) => source.len_hint(),
        }
    }
}

#[cfg(feature = "gzip")]
fn open_gzip(
    path: &Path,
    window_bytes: usize,
    context: &ExecutionContext,
) -> Result<SpillWindowedFile, Diagnostic> {
    let file = File::open(path).map_err(|error| io_error("open gzip", &error))?;
    spool(flate2::read::GzDecoder::new(file), window_bytes, context)
}

#[cfg(not(feature = "gzip"))]
fn open_gzip(
    _path: &Path,
    _window_bytes: usize,
    _context: &ExecutionContext,
) -> Result<SpillWindowedFile, Diagnostic> {
    Err(codec_error("gzip support is not enabled"))
}

#[cfg(feature = "zstd")]
fn open_zstd(
    path: &Path,
    window_bytes: usize,
    context: &ExecutionContext,
) -> Result<SpillWindowedFile, Diagnostic> {
    let file = File::open(path).map_err(|error| io_error("open zstd", &error))?;
    let decoder = zstd::stream::read::Decoder::new(file)
        .map_err(|error| codec_error(format!("zstd decoder failed: {error}")))?;
    spool(decoder, window_bytes, context)
}

#[cfg(not(feature = "zstd"))]
fn open_zstd(
    _path: &Path,
    _window_bytes: usize,
    _context: &ExecutionContext,
) -> Result<SpillWindowedFile, Diagnostic> {
    Err(codec_error("zstd support is not enabled"))
}

#[cfg(any(feature = "gzip", feature = "zstd"))]
fn spool(
    mut reader: impl Read,
    record_bytes: usize,
    context: &ExecutionContext,
) -> Result<SpillWindowedFile, Diagnostic> {
    let key = format!(
        "structure-input-{}",
        SPOOL_SEQUENCE.fetch_add(1, Ordering::Relaxed)
    );
    let mut spill = context
        .create_spill_file(&key)
        .map_err(|error| spill_error(&error))?;
    let reservation = context
        .try_reserve(record_bytes)
        .map_err(|error| resource_error(error.to_string()))?;
    let mut buffer = Vec::new();
    buffer
        .try_reserve_exact(record_bytes)
        .map_err(|error| resource_error(error.to_string()))?;
    buffer.resize(record_bytes, 0);
    let mut length = 0_u64;
    loop {
        let mut filled = 0usize;
        while filled < record_bytes {
            match reader.read(&mut buffer[filled..]) {
                Ok(0) => break,
                Ok(count) => filled += count,
                Err(error) if error.kind() == ErrorKind::Interrupted => {}
                Err(error) => return Err(codec_error(format!("decompression failed: {error}"))),
            }
        }
        if filled == 0 {
            break;
        }
        spill
            .append(&buffer[..filled])
            .map_err(|error| spill_error(&error))?;
        length = length
            .checked_add(
                u64::try_from(filled)
                    .map_err(|_| resource_error("decompressed length exceeds u64"))?,
            )
            .ok_or_else(|| resource_error("decompressed length overflow"))?;
        if filled < record_bytes {
            break;
        }
    }
    let artifact = spill.finish().map_err(|error| spill_error(&error))?;
    drop(buffer);
    drop(reservation);
    SpillWindowedFile::open(artifact, record_bytes, length, record_bytes, context)
}

fn io_error(operation: &'static str, error: &std::io::Error) -> Diagnostic {
    Diagnostic::new(Code::E1901)
        .with_message("structural source I/O failed")
        .with_context("operation", operation)
        .with_context("reason", error.to_string())
}

fn codec_error(reason: impl Into<Box<str>>) -> Diagnostic {
    Diagnostic::new(Code::E1401).with_message(reason)
}

#[cfg(any(feature = "gzip", feature = "zstd"))]
fn spill_error(error: &crate::SpillError) -> Diagnostic {
    resource_error(format!(
        "compressed structure requires bounded spill: {error}"
    ))
}

#[cfg(any(feature = "gzip", feature = "zstd"))]
fn resource_error(reason: impl Into<Box<str>>) -> Diagnostic {
    Diagnostic::new(Code::E1902).with_message(reason)
}

#[cfg(test)]
#[path = "source_file_tests.rs"]
mod tests;
