//! Direct frame positioning for file-backed DCD sources.
//!
//! Every frame occupies a computable number of bytes: the first carries all
//! atoms and later ones only the free atoms when some are fixed. Seeking sums
//! those sizes instead of decoding, so random access costs one file seek.

use std::io::{Seek, SeekFrom};
use std::mem::size_of;

use super::{DcdReader, TrajectoryError, io_error, trajectory_error};
use crate::Timestep;

const CELL_RECORD_BYTES: u64 = 48 + 8;
const MARKER_BYTES: u64 = 8;

impl DcdReader {
    pub(super) fn seek_frame(&mut self, frame: usize) -> Result<(), TrajectoryError> {
        if frame > self.header.frame_count {
            return Err(TrajectoryError::InvalidSource {
                format: "dcd frame out of range",
            });
        }
        self.load_fixed_atoms().map_err(|e| trajectory_error(&e))?;
        let position = self.frame_position(frame)?;
        self.records
            .input
            .seek(SeekFrom::Start(position))
            .map_err(|error| trajectory_error(&io_error(&error)))?;
        self.records.pending_length = None;
        self.records.offset =
            usize::try_from(position).map_err(|_| TrajectoryError::InvalidSource {
                format: "dcd frame out of range",
            })?;
        self.next_frame = frame;
        self.eof_checked = false;
        Ok(())
    }

    /// Fixed atoms are stored only in frame zero, so later frames need it first.
    fn load_fixed_atoms(&mut self) -> Result<(), super::DcdError> {
        if self.header.fixed_atom_count == 0
            || self.fixed_positions.len() == self.header.fixed_atom_count
            || self.header.frame_count == 0
        {
            return Ok(());
        }
        self.records
            .input
            .seek(SeekFrom::Start(self.data_start))
            .map_err(|error| io_error(&error))?;
        self.records.pending_length = None;
        self.next_frame = 0;
        let mut scratch = Timestep::default();
        self.read_next_frame(&mut scratch)?;
        Ok(())
    }

    fn frame_position(&self, frame: usize) -> Result<u64, TrajectoryError> {
        let overflow = TrajectoryError::InvalidSource {
            format: "dcd frame out of range",
        };
        let first = self.frame_bytes(self.header.atom_count);
        let later = self.frame_bytes(self.header.atom_count - self.header.fixed_atom_count);
        let skipped = match frame {
            0 => Some(0),
            _ => u64::try_from(frame - 1)
                .ok()
                .and_then(|count| count.checked_mul(later))
                .and_then(|bytes| bytes.checked_add(first)),
        };
        skipped
            .and_then(|bytes| bytes.checked_add(self.data_start))
            .ok_or(overflow)
    }

    fn frame_bytes(&self, atoms: usize) -> u64 {
        let axis = (atoms as u64) * (size_of::<f32>() as u64) + MARKER_BYTES;
        let cell = if self.header.has_unit_cell {
            CELL_RECORD_BYTES
        } else {
            0
        };
        let fourth = if self.header.has_fourth_dimension {
            axis
        } else {
            0
        };
        cell + 3 * axis + fourth
    }
}

#[cfg(test)]
#[path = "reader_seek_tests.rs"]
mod tests;
