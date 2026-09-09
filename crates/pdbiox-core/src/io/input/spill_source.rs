//! Random-access windows over fixed-record execution-owned spill storage.

use super::{ByteWindow, SourceBytes};
use crate::execution::spill::{FILE_HEADER_BYTES, RECORD_HEADER_DISK_BYTES, checksum};
use crate::{Code, Diagnostic, ExecutionContext, SpillArtifact};
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};

const HEADER_BYTES: usize = 16;

/// A bounded source backed by an immutable, automatically deleted spill file.
#[derive(Debug)]
pub struct SpillWindowedFile {
    file: File,
    artifact: SpillArtifact,
    record_bytes: usize,
    length: u64,
    output: Vec<u8>,
    record: Vec<u8>,
    cached_record: Option<u64>,
    _reservation: crate::MemoryReservation,
}

impl SpillWindowedFile {
    /// Opens fixed-size records as one logical byte stream.
    ///
    /// # Errors
    ///
    /// Returns a resource diagnostic if buffers, the artifact, or its record
    /// geometry are invalid.
    pub fn open(
        artifact: SpillArtifact,
        record_bytes: usize,
        length: u64,
        max_window_bytes: usize,
        context: &ExecutionContext,
    ) -> Result<Self, Diagnostic> {
        if record_bytes == 0 || max_window_bytes == 0 {
            return Err(resource_error("spill source buffers must be non-zero"));
        }
        let retained = record_bytes
            .checked_add(max_window_bytes)
            .ok_or_else(|| resource_error("spill source buffer size overflow"))?;
        let reservation = context
            .try_reserve(retained)
            .map_err(|error| resource_error(error.to_string()))?;
        let file = File::open(artifact.path()).map_err(|error| io_error("open", &error))?;
        let mut output = Vec::new();
        output
            .try_reserve_exact(max_window_bytes)
            .map_err(|error| resource_error(error.to_string()))?;
        let mut record = Vec::new();
        record
            .try_reserve_exact(record_bytes)
            .map_err(|error| resource_error(error.to_string()))?;
        Ok(Self {
            file,
            artifact,
            record_bytes,
            length,
            output,
            record,
            cached_record: None,
            _reservation: reservation,
        })
    }

    fn load_record(&mut self, index: u64) -> Result<(), Diagnostic> {
        if self.cached_record == Some(index) {
            return Ok(());
        }
        let stride = u64::try_from(self.record_bytes)
            .map_err(|_| resource_error("spill record size exceeds u64"))?
            .checked_add(RECORD_HEADER_DISK_BYTES)
            .ok_or_else(|| resource_error("spill record stride overflow"))?;
        let offset = index
            .checked_mul(stride)
            .and_then(|value| value.checked_add(FILE_HEADER_BYTES))
            .ok_or_else(|| resource_error("spill record offset overflow"))?;
        self.file
            .seek(SeekFrom::Start(offset))
            .map_err(|error| io_error("seek", &error))?;
        let mut header = [0_u8; HEADER_BYTES];
        self.file
            .read_exact(&mut header)
            .map_err(|error| io_error("read header", &error))?;
        let length = u64::from_le_bytes(copy_u64(&header[..8]));
        let expected = u64::from_le_bytes(copy_u64(&header[8..]));
        let maximum = u64::try_from(self.record_bytes)
            .map_err(|_| resource_error("spill record size exceeds u64"))?;
        if length > maximum {
            return Err(resource_error("spill record exceeds configured width"));
        }
        let length =
            usize::try_from(length).map_err(|_| resource_error("spill record exceeds usize"))?;
        self.record.resize(length, 0);
        self.file
            .read_exact(&mut self.record)
            .map_err(|error| io_error("read payload", &error))?;
        if checksum(&self.record) != expected {
            return Err(resource_error("spill record checksum mismatch"));
        }
        self.cached_record = Some(index);
        Ok(())
    }
}

impl SourceBytes for SpillWindowedFile {
    fn window(&mut self, start: u64, len: usize) -> Result<ByteWindow<'_>, Diagnostic> {
        if len > self.output.capacity() {
            return Err(resource_error("requested spill window exceeds capacity"));
        }
        if start >= self.length || len == 0 {
            return Ok(ByteWindow::new(start, &[]));
        }
        let remaining = usize::try_from(self.length - start).map_or(usize::MAX, |value| value);
        let target = len.min(remaining);
        self.output.clear();
        let record_bytes = u64::try_from(self.record_bytes)
            .map_err(|_| resource_error("spill record size exceeds u64"))?;
        let mut position = start;
        while self.output.len() < target {
            let index = position / record_bytes;
            let local = usize::try_from(position % record_bytes)
                .map_err(|_| resource_error("spill record offset exceeds usize"))?;
            self.load_record(index)?;
            let available = self.record.len().saturating_sub(local);
            if available == 0 {
                return Err(resource_error(
                    "spill source ended before its logical length",
                ));
            }
            let count = available.min(target.saturating_sub(self.output.len()));
            let end = local.saturating_add(count);
            let bytes = self
                .record
                .get(local..end)
                .ok_or_else(|| resource_error("spill record range is invalid"))?;
            self.output.extend_from_slice(bytes);
            position = position
                .checked_add(
                    u64::try_from(count).map_err(|_| resource_error("spill copy exceeds u64"))?,
                )
                .ok_or_else(|| resource_error("spill source offset overflow"))?;
        }
        Ok(ByteWindow::new(start, &self.output))
    }

    fn len_hint(&self) -> Option<u64> {
        Some(self.length)
    }
}

impl SpillWindowedFile {
    /// Physical bytes retained by the execution-owned artifact.
    #[must_use]
    pub const fn spill_bytes(&self) -> u64 {
        self.artifact.bytes()
    }
}

fn copy_u64(bytes: &[u8]) -> [u8; 8] {
    let mut value = [0_u8; 8];
    value.copy_from_slice(bytes);
    value
}

fn io_error(operation: &'static str, error: &std::io::Error) -> Diagnostic {
    resource_error(format!("spill source {operation} failed: {error}"))
}

fn resource_error(reason: impl Into<Box<str>>) -> Diagnostic {
    Diagnostic::new(Code::E1902).with_message(reason)
}
