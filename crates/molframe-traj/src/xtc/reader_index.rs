//! Sharing and persisting the frame-offset index of an XTC reader.

use std::path::Path;

use super::{XtcError, XtcReader, memory_error};
use crate::frame_index::{FrameIndex, FrameIndexError};

impl XtcReader {
    /// The frame-offset index, built by one header scan if not yet known.
    ///
    /// The scan reads only frame headers and the payload length word, so its
    /// cost is bounded by the number of frames, not the file size.
    ///
    /// # Errors
    ///
    /// Returns an I/O, malformed-record or memory-ceiling error.
    pub fn frame_index(&mut self) -> Result<FrameIndex, XtcError> {
        self.ensure_index()?;
        let offsets = self.offsets.as_ref().ok_or(XtcError::InvalidFrame)?;
        FrameIndex::new(
            FrameIndex::XTC,
            self.atom_count,
            self.stamp,
            offsets.to_vec(),
        )
        .map_err(XtcError::from)
    }

    /// Adopts an index built elsewhere, skipping this reader's own scan.
    ///
    /// # Errors
    ///
    /// Returns [`FrameIndexError::Mismatch`] (as [`XtcError::Index`]) when the
    /// index describes another file, format or atom count, or a memory error
    /// when the index does not fit the ceiling.
    pub fn with_frame_index(mut self, index: &FrameIndex) -> Result<Self, XtcError> {
        self.adopt_frame_index(index)?;
        Ok(self)
    }

    /// In-place form of [`Self::with_frame_index`].
    ///
    /// # Errors
    ///
    /// As [`Self::with_frame_index`].
    pub fn adopt_frame_index(&mut self, index: &FrameIndex) -> Result<(), XtcError> {
        if !index.describes(FrameIndex::XTC, self.atom_count, self.stamp) {
            return Err(FrameIndexError::Mismatch.into());
        }
        let required = self.total_bytes(self.atom_count, self.scratch.capacity(), index.len());
        if required > self.memory_limit_bytes {
            return Err(memory_error(required, self.memory_limit_bytes));
        }
        self.frame_count = Some(index.len());
        self.offsets = Some(index.shared());
        Ok(())
    }

    /// Opens a file with its index loaded from `cache` when that still matches.
    ///
    /// A missing, damaged or stale cache is rebuilt by a header scan and
    /// rewritten; a cache that cannot be written is ignored, because it only
    /// ever saves time.
    ///
    /// # Errors
    ///
    /// Returns the open and scan errors of [`Self::open_with_memory_limit`]
    /// and [`Self::frame_index`].
    pub fn open_cached(
        path: &Path,
        memory_limit_bytes: usize,
        cache: &Path,
    ) -> Result<Self, XtcError> {
        let reader = Self::open_with_memory_limit(path, memory_limit_bytes)?;
        if let Ok(index) = FrameIndex::read(cache)
            && index.describes(FrameIndex::XTC, reader.atom_count, reader.stamp)
        {
            return reader.with_frame_index(&index);
        }
        let mut reader = reader;
        let index = reader.frame_index()?;
        let _ = index.write(cache);
        Ok(reader)
    }
}

#[cfg(test)]
#[path = "reader_index_tests.rs"]
mod tests;
