//! Frame-offset index, seek and index sharing for the TRR reader.

use std::fs::File;
use std::io::{BufReader, Seek, SeekFrom};
use std::sync::Arc;

use super::{BUFFER, TrrReader, invalid, io_error, memory, next_header};
use crate::frame_index::FrameIndex;
use crate::{FileStamp, TrajectoryError};

impl TrrReader {
    /// The frame-offset index, built by one header scan if not yet known.
    ///
    /// The scan runs on its own file handle and reads only frame headers, so
    /// it does not disturb the reader's position.
    ///
    /// # Errors
    ///
    /// Returns a malformed-header, I/O or memory error.
    pub fn frame_index(&mut self) -> Result<FrameIndex, TrajectoryError> {
        self.ensure_index()?;
        let offsets = self.offsets.as_ref().ok_or_else(invalid)?;
        Ok(FrameIndex::new(
            FrameIndex::TRR,
            self.atoms,
            self.stamp,
            offsets.to_vec(),
        )?)
    }

    /// Adopts an index built elsewhere, skipping this reader's own scan.
    ///
    /// # Errors
    ///
    /// Returns an invalid-source error when the index describes another file,
    /// format or atom count, or a memory error when it exceeds the ceiling.
    pub fn with_frame_index(mut self, index: &FrameIndex) -> Result<Self, TrajectoryError> {
        self.adopt_frame_index(index)?;
        Ok(self)
    }

    /// In-place form of [`Self::with_frame_index`].
    ///
    /// # Errors
    ///
    /// As [`Self::with_frame_index`].
    pub fn adopt_frame_index(&mut self, index: &FrameIndex) -> Result<(), TrajectoryError> {
        if !index.describes(FrameIndex::TRR, self.atoms, self.stamp) {
            return Err(crate::FrameIndexError::Mismatch.into());
        }
        let required = index.len().saturating_mul(size_of::<u64>());
        if required.saturating_add(BUFFER) > self.limit {
            return Err(memory(required, self.limit));
        }
        self.offsets = Some(index.shared());
        Ok(())
    }

    pub(super) fn seek_indexed(&mut self, frame: usize) -> Result<(), TrajectoryError> {
        self.ensure_index()?;
        let offsets = self.offsets.as_ref().ok_or_else(invalid)?;
        let target = frame.min(offsets.len());
        let offset = offsets.get(target).map_or(self.file_len, |at| *at);
        let previous = target
            .checked_sub(1)
            .and_then(|before| offsets.get(before))
            .copied();
        // The step that precedes the target supplies `dt`, so a frame decoded
        // after a seek carries what a sequential read would have given it.
        let previous_time = match previous {
            Some(before) => {
                self.input
                    .seek(SeekFrom::Start(before))
                    .map_err(|error| io_error(&error))?;
                next_header(&mut self.input)?.map(|header| header.time)
            }
            None => None,
        };
        self.input
            .seek(SeekFrom::Start(offset))
            .map_err(|error| io_error(&error))?;
        self.pending = None;
        self.next = target;
        self.previous_time = previous_time;
        self.last_step = None;
        self.last_precision = None;
        Ok(())
    }

    fn ensure_index(&mut self) -> Result<(), TrajectoryError> {
        if self.offsets.is_some() {
            return Ok(());
        }
        let file = File::open(&self.path).map_err(|error| io_error(&error))?;
        let stamp = FileStamp::from_metadata(&file.metadata().map_err(|error| io_error(&error))?);
        if stamp != self.stamp {
            return Err(invalid());
        }
        let mut input = BufReader::with_capacity(BUFFER, file);
        let mut offsets = Vec::new();
        let mut cursor = 0_u64;
        while cursor < self.file_len {
            let header = next_header(&mut input)?.ok_or_else(invalid)?;
            if header.atoms != self.atoms {
                return Err(invalid());
            }
            let payload = header
                .sizes
                .iter()
                .try_fold(0_u64, |sum, size| sum.checked_add(*size as u64))
                .ok_or_else(invalid)?;
            let after = input.stream_position().map_err(|error| io_error(&error))?;
            let end = after.checked_add(payload).ok_or_else(invalid)?;
            if end > self.file_len {
                return Err(invalid());
            }
            self.reserve_slot(&mut offsets)?;
            offsets.push(cursor);
            input
                .seek(SeekFrom::Start(end))
                .map_err(|error| io_error(&error))?;
            cursor = end;
        }
        self.offsets = Some(Arc::from(offsets));
        Ok(())
    }

    fn reserve_slot(&self, offsets: &mut Vec<u64>) -> Result<(), TrajectoryError> {
        let needed = (offsets.len() + 1).saturating_mul(size_of::<u64>());
        if needed.saturating_add(BUFFER) > self.limit {
            return Err(memory(needed, self.limit));
        }
        offsets
            .try_reserve(1)
            .map_err(|_| memory(usize::MAX, self.limit))
    }
}
