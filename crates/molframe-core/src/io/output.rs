//! Deterministic, bounded output selected by the destination suffix.

use crate::diagnostic::{Code, Diagnostic};
use std::fmt;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

/// Default ceiling for format-owned output workspace.
///
/// A default, not a maximum. A writer's workspace is bounded by what the caller
/// allows, and a caller writing a very large structure may allow more.
pub const DEFAULT_OUTPUT_MEMORY_LIMIT_BYTES: usize = 100_000_000;
const OUTPUT_BUFFER_BYTES: usize = 64 * 1_024;

/// Working-memory policy for incremental output.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct OutputOptions {
    /// Maximum bytes a format writer may retain as output workspace.
    pub memory_limit_bytes: usize,
}

impl OutputOptions {
    /// Replaces the default 100 MB output-workspace ceiling.
    #[must_use]
    pub const fn with_memory_limit(mut self, bytes: usize) -> Self {
        self.memory_limit_bytes = bytes;
        self
    }
}

impl Default for OutputOptions {
    fn default() -> Self {
        Self {
            memory_limit_bytes: DEFAULT_OUTPUT_MEMORY_LIMIT_BYTES,
        }
    }
}

/// Adapts a byte writer for allocation-free [`std::fmt::Write`] formatting.
///
/// Formatting APIs expose only `fmt::Error`; this adapter retains the original
/// I/O error so callers can recover it with [`TextOutput::finish`].
pub struct TextOutput<'a, W> {
    output: &'a mut W,
    error: Option<io::Error>,
}

impl<W> fmt::Debug for TextOutput<'_, W> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("TextOutput")
            .field("failed", &self.error.is_some())
            .finish_non_exhaustive()
    }
}

impl<'a, W> TextOutput<'a, W> {
    /// Borrows a byte destination without adding a buffer.
    pub const fn new(output: &'a mut W) -> Self {
        Self {
            output,
            error: None,
        }
    }

    /// Returns the first destination failure, if any.
    ///
    /// # Errors
    ///
    /// Returns the original I/O error retained while formatting.
    pub fn finish(self) -> io::Result<()> {
        match self.error {
            Some(error) => Err(error),
            None => Ok(()),
        }
    }
}

impl<W: Write> std::fmt::Write for TextOutput<'_, W> {
    fn write_str(&mut self, text: &str) -> std::fmt::Result {
        if self.error.is_some() {
            return Err(std::fmt::Error);
        }
        if let Err(error) = self.output.write_all(text.as_bytes()) {
            self.error = Some(error);
            return Err(std::fmt::Error);
        }
        Ok(())
    }
}

/// A deterministic compressed or plain output stream.
///
/// The sink owns at most 64 KiB of Rust buffering. Format writers must keep
/// their additional workspace within [`OutputOptions::memory_limit_bytes`].
#[must_use = "compressed streams must be finished explicitly"]
pub struct OutputSink<W: Write> {
    encoder: Option<Encoder<W>>,
    destination: Box<str>,
}

impl<W: Write> fmt::Debug for OutputSink<W> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("OutputSink")
            .field("destination", &self.destination)
            .field("finished", &self.encoder.is_none())
            .finish_non_exhaustive()
    }
}

enum Encoder<W: Write> {
    Plain(BufWriter<W>),
    #[cfg(feature = "gzip")]
    Gzip(flate2::write::GzEncoder<BufWriter<W>>),
    #[cfg(feature = "zstd")]
    Zstd(zstd::stream::write::Encoder<'static, BufWriter<W>>),
}

impl OutputSink<File> {
    /// Creates a local output stream, selecting compression from its suffix.
    ///
    /// # Errors
    ///
    /// Returns `E7901` for invalid limits, unavailable compression, or an
    /// output file that cannot be created.
    pub fn create(path: impl AsRef<Path>, options: OutputOptions) -> Result<Self, Diagnostic> {
        let path = path.as_ref();
        validate_options(path, options)?;
        let compression = compression_for(path);
        ensure_available(path, compression)?;
        let file = File::create(path).map_err(|error| output_error(path, error))?;
        Self::from_writer(
            file,
            compression,
            options,
            path.display().to_string().into(),
        )
        .map_err(|error| output_error(path, error))
    }
}

impl<W: Write> OutputSink<W> {
    fn from_writer(
        writer: W,
        compression: OutputCompression,
        options: OutputOptions,
        destination: Box<str>,
    ) -> io::Result<Self> {
        let capacity = OUTPUT_BUFFER_BYTES.min(options.memory_limit_bytes);
        let writer = BufWriter::with_capacity(capacity, writer);
        let encoder = match compression {
            OutputCompression::None => Encoder::Plain(writer),
            #[cfg(feature = "gzip")]
            OutputCompression::Gzip => Encoder::Gzip(
                flate2::GzBuilder::new()
                    .mtime(0)
                    .write(writer, flate2::Compression::default()),
            ),
            #[cfg(feature = "zstd")]
            OutputCompression::Zstd => Encoder::Zstd(zstd::stream::write::Encoder::new(writer, 0)?),
            #[cfg(any(not(feature = "gzip"), not(feature = "zstd")))]
            _ => return Err(io::Error::other("compression support is not enabled")),
        };
        Ok(Self {
            encoder: Some(encoder),
            destination,
        })
    }

    /// Completes compression and flushes all bytes to the underlying writer.
    ///
    /// # Errors
    ///
    /// Returns `E7901` if the stream cannot be completed.
    pub fn finish(mut self) -> Result<(), Diagnostic> {
        let Some(encoder) = self.encoder.take() else {
            return Err(output_error_name(
                &self.destination,
                "output stream was already finished",
            ));
        };
        finish_encoder(encoder).map_err(|error| output_error_name(&self.destination, error))
    }

    #[cfg(test)]
    fn buffer_capacity(&self) -> Option<usize> {
        match self.encoder.as_ref() {
            Some(Encoder::Plain(writer)) => Some(writer.capacity()),
            #[cfg(feature = "gzip")]
            Some(Encoder::Gzip(writer)) => Some(writer.get_ref().capacity()),
            #[cfg(feature = "zstd")]
            Some(Encoder::Zstd(writer)) => Some(writer.get_ref().capacity()),
            None => None,
        }
    }
}

impl<W: Write> Write for OutputSink<W> {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        match self.encoder.as_mut() {
            Some(Encoder::Plain(writer)) => writer.write(bytes),
            #[cfg(feature = "gzip")]
            Some(Encoder::Gzip(writer)) => writer.write(bytes),
            #[cfg(feature = "zstd")]
            Some(Encoder::Zstd(writer)) => writer.write(bytes),
            None => Err(io::Error::new(
                io::ErrorKind::BrokenPipe,
                "output stream is finished",
            )),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self.encoder.as_mut() {
            Some(Encoder::Plain(writer)) => writer.flush(),
            #[cfg(feature = "gzip")]
            Some(Encoder::Gzip(writer)) => writer.flush(),
            #[cfg(feature = "zstd")]
            Some(Encoder::Zstd(writer)) => writer.flush(),
            None => Err(io::Error::new(
                io::ErrorKind::BrokenPipe,
                "output stream is finished",
            )),
        }
    }
}

/// Writes existing bytes through the bounded incremental sink.
///
/// Prefer [`OutputSink`] when the producer can generate its output directly;
/// this compatibility helper does not allocate another full-file buffer.
///
/// # Errors
///
/// Returns `E7901` when compression is unavailable or output fails.
pub fn write_output(path: impl AsRef<Path>, bytes: &[u8]) -> Result<(), Diagnostic> {
    let path = path.as_ref();
    let mut output = OutputSink::create(path, OutputOptions::default())?;
    output
        .write_all(bytes)
        .map_err(|error| output_error(path, error))?;
    output.finish()
}

#[cfg(any(feature = "gzip", feature = "zstd"))]
fn finish_encoder<W: Write>(encoder: Encoder<W>) -> io::Result<()> {
    let mut writer = match encoder {
        Encoder::Plain(writer) => writer,
        #[cfg(feature = "gzip")]
        Encoder::Gzip(writer) => writer.finish()?,
        #[cfg(feature = "zstd")]
        Encoder::Zstd(writer) => writer.finish()?,
    };
    writer.flush()
}

#[cfg(not(any(feature = "gzip", feature = "zstd")))]
fn finish_encoder<W: Write>(encoder: Encoder<W>) -> io::Result<()> {
    let Encoder::Plain(mut writer) = encoder;
    writer.flush()
}

#[derive(Clone, Copy)]
enum OutputCompression {
    None,
    Gzip,
    Zstd,
}

fn compression_for(path: &Path) -> OutputCompression {
    match path.extension().and_then(|suffix| suffix.to_str()) {
        Some(suffix) if suffix.eq_ignore_ascii_case("gz") => OutputCompression::Gzip,
        Some(suffix) if suffix.eq_ignore_ascii_case("zst") => OutputCompression::Zstd,
        _ => OutputCompression::None,
    }
}

// Every arm is compiled out when both codecs are enabled, leaving a function
// that cannot fail. The signature stays uniform across feature sets so the call
// site does not change shape with the build.
#[cfg_attr(
    all(feature = "gzip", feature = "zstd"),
    allow(clippy::unnecessary_wraps)
)]
fn ensure_available(path: &Path, compression: OutputCompression) -> Result<(), Diagnostic> {
    let _ = path;
    match compression {
        #[cfg(not(feature = "gzip"))]
        OutputCompression::Gzip => Err(unavailable(path, "gzip")),
        #[cfg(not(feature = "zstd"))]
        OutputCompression::Zstd => Err(unavailable(path, "zstandard")),
        _ => Ok(()),
    }
}

fn validate_options(path: &Path, options: OutputOptions) -> Result<(), Diagnostic> {
    if options.memory_limit_bytes == 0 {
        return Err(Diagnostic::new(Code::E7901)
            .with_context("path", path.display().to_string())
            .with_context("memory limit", options.memory_limit_bytes.to_string()));
    }
    Ok(())
}

#[cfg(any(not(feature = "gzip"), not(feature = "zstd")))]
fn unavailable(path: &Path, compression: &str) -> Diagnostic {
    Diagnostic::new(Code::E7901)
        .with_context("path", path.display().to_string())
        .with_context("compression", compression)
        .with_context("reason", "support is not enabled in this build")
}

fn output_error(path: &Path, error: impl std::fmt::Display) -> Diagnostic {
    output_error_name(&path.display().to_string(), error)
}

fn output_error_name(destination: &str, error: impl std::fmt::Display) -> Diagnostic {
    Diagnostic::new(Code::E7901)
        .with_context("path", destination)
        .with_context("reason", error.to_string())
}

#[cfg(test)]
#[path = "output_tests.rs"]
mod tests;
