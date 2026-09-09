//! Stable bounded windows over resident and file-backed sources.

use super::InputBuffer;
use crate::{Code, Diagnostic, ExecutionContext, MemoryReservation};
use std::fs::File;
use std::io::{ErrorKind, Read, Seek, SeekFrom};
use std::path::Path;

/// One stable view into a source.
#[derive(Clone, Copy, Debug)]
pub struct ByteWindow<'a> {
    start: u64,
    bytes: &'a [u8],
}

impl<'a> ByteWindow<'a> {
    /// Creates a view whose first byte has the supplied global offset.
    #[must_use]
    pub const fn new(start: u64, bytes: &'a [u8]) -> Self {
        Self { start, bytes }
    }

    /// Global offset of the first byte.
    #[must_use]
    pub const fn start(self) -> u64 {
        self.start
    }

    /// Bytes present in this window. A final window may be shorter than demand.
    #[must_use]
    pub const fn bytes(self) -> &'a [u8] {
        self.bytes
    }

    /// Exclusive global end offset.
    ///
    /// # Errors
    ///
    /// Returns `PDBIOX-E1903` if the exclusive end would exceed `u64::MAX`.
    pub fn end(self) -> Result<u64, Diagnostic> {
        let length =
            u64::try_from(self.bytes.len()).map_err(|_| unrepresentable(self.bytes.len()))?;
        self.start
            .checked_add(length)
            .ok_or_else(|| unrepresentable("an offset greater than u64::MAX"))
    }
}

/// Random-access bytes exposed through bounded, borrow-stable windows.
pub trait SourceBytes {
    /// Returns at most `len` bytes starting at the global byte offset `start`.
    ///
    /// # Errors
    ///
    /// Returns a registered resource, budget, or addressability diagnostic.
    fn window(&mut self, start: u64, len: usize) -> Result<ByteWindow<'_>, Diagnostic>;

    /// Exact source length when known without consuming the source.
    fn len_hint(&self) -> Option<u64>;
}

impl SourceBytes for InputBuffer {
    fn window(&mut self, start: u64, len: usize) -> Result<ByteWindow<'_>, Diagnostic> {
        let local = usize::try_from(start).map_err(|_| unrepresentable(start))?;
        let bytes = self.as_bytes();
        if local >= bytes.len() {
            return Ok(ByteWindow::new(start, &[]));
        }
        let end = local.saturating_add(len).min(bytes.len());
        Ok(ByteWindow::new(start, &bytes[local..end]))
    }

    fn len_hint(&self) -> Option<u64> {
        u64::try_from(self.len()).ok()
    }
}

/// A seekable file read through one reusable, budget-charged window.
#[derive(Debug)]
pub struct WindowedFile {
    file: File,
    length: u64,
    buffer: Box<[u8]>,
    _reservation: MemoryReservation,
}

impl WindowedFile {
    /// Opens a file and reserves the complete reusable window up front.
    ///
    /// # Errors
    ///
    /// Returns `PDBIOX-E1902` if the window does not fit the shared execution
    /// budget, or `PDBIOX-E1901` if the file or buffer cannot be opened.
    pub fn open(
        path: impl AsRef<Path>,
        max_window_bytes: usize,
        context: &ExecutionContext,
    ) -> Result<Self, Diagnostic> {
        if max_window_bytes == 0 {
            return Err(Diagnostic::new(Code::E1902).with_context("requested bytes", "0"));
        }
        let reservation = context.try_reserve(max_window_bytes).map_err(|error| {
            Diagnostic::new(Code::E1902).with_context("reason", error.to_string())
        })?;
        let file = File::open(path.as_ref())
            .map_err(|error| io_failure("source file could not be opened", &error))?;
        let length = file
            .metadata()
            .map_err(|error| io_failure("source length could not be read", &error))?
            .len();
        let mut buffer = Vec::new();
        buffer
            .try_reserve_exact(max_window_bytes)
            .map_err(|error| {
                Diagnostic::new(Code::E1901).with_context("reason", error.to_string())
            })?;
        buffer.resize(max_window_bytes, 0);
        Ok(Self {
            file,
            length,
            buffer: buffer.into_boxed_slice(),
            _reservation: reservation,
        })
    }

    /// Maximum resident source-window capacity.
    #[must_use]
    pub const fn capacity(&self) -> usize {
        self.buffer.len()
    }
}

impl SourceBytes for WindowedFile {
    fn window(&mut self, start: u64, len: usize) -> Result<ByteWindow<'_>, Diagnostic> {
        if len > self.buffer.len() {
            return Err(Diagnostic::new(Code::E1902)
                .with_context("requested bytes", len.to_string())
                .with_context("window bytes", self.buffer.len().to_string()));
        }
        if start >= self.length || len == 0 {
            return Ok(ByteWindow::new(start, &[]));
        }
        self.file
            .seek(SeekFrom::Start(start))
            .map_err(|error| io_failure("source window could not seek", &error))?;
        let remaining = self.length - start;
        let host_remaining = match usize::try_from(remaining) {
            Ok(value) => value,
            Err(_) => usize::MAX,
        };
        let target = len.min(host_remaining);
        let mut filled = 0;
        while filled < target {
            match self.file.read(&mut self.buffer[filled..target]) {
                Ok(0) => break,
                Ok(count) => filled += count,
                Err(error) if error.kind() == ErrorKind::Interrupted => {}
                Err(error) => return Err(io_failure("source window could not be read", &error)),
            }
        }
        Ok(ByteWindow::new(start, &self.buffer[..filled]))
    }

    fn len_hint(&self) -> Option<u64> {
        Some(self.length)
    }
}

fn io_failure(message: &'static str, error: &std::io::Error) -> Diagnostic {
    Diagnostic::new(Code::E1901)
        .with_message(message)
        .with_context("reason", error.to_string())
}

fn unrepresentable(value: impl std::fmt::Display) -> Diagnostic {
    Diagnostic::new(Code::E1903).with_context("value", value.to_string())
}

#[cfg(test)]
#[path = "source_tests.rs"]
mod tests;
