//! Common trajectory-reader operations for the bounded XTC decoder.

use super::{
    FORMAT, XtcReader,
    bounded::{guarded_record, trajectory_error},
};
use crate::{RandomAccess, Timestep, TrajectoryError, TrajectoryReader, Units};
use std::io::{Seek, SeekFrom};

impl TrajectoryReader for XtcReader {
    fn format(&self) -> &'static str {
        FORMAT
    }

    fn n_atoms(&self) -> usize {
        self.atom_count
    }

    fn n_frames(&self) -> Option<usize> {
        self.frame_count
    }

    fn units(&self) -> Units {
        Units::CANONICAL
    }

    fn random_access(&self) -> RandomAccess {
        if self.offsets.is_some() {
            RandomAccess::Full
        } else if self.index_unavailable {
            RandomAccess::None
        } else {
            RandomAccess::ViaIndex
        }
    }

    fn read_next(&mut self, timestep: &mut Timestep) -> Result<bool, TrajectoryError> {
        self.read_next_frame(timestep).map_err(trajectory_error)
    }

    fn read_next_bounded(
        &mut self,
        timestep: &mut Timestep,
        bytes: usize,
    ) -> Result<bool, TrajectoryError> {
        let previous = self.memory_limit_bytes;
        self.memory_limit_bytes = previous.min(bytes);
        let result = self.read_next(timestep);
        self.memory_limit_bytes = previous;
        result
    }

    fn seek(&mut self, frame: usize) -> Result<(), TrajectoryError> {
        if self.index_unavailable {
            return Err(TrajectoryError::RandomAccessUnavailable);
        }
        self.ensure_index().map_err(trajectory_error)?;
        let offsets = self
            .offsets
            .as_ref()
            .ok_or(TrajectoryError::RandomAccessUnavailable)?;
        let target = frame.min(offsets.len());
        let offset = match offsets.get(target) {
            Some(offset) => *offset,
            None => self.file_len,
        };
        let before = target
            .checked_sub(1)
            .and_then(|previous| offsets.get(previous))
            .copied();
        let io = |error: std::io::Error| TrajectoryError::SourceIo {
            format: FORMAT,
            kind: error.kind(),
        };
        // The header before the target supplies `dt`, so a frame decoded after
        // a seek carries what a sequential read would have given it.
        let previous_time = match before {
            Some(before) => {
                self.source.file.seek(SeekFrom::Start(before)).map_err(io)?;
                let header =
                    guarded_record(|| self.source.read_header()).map_err(trajectory_error)?;
                Some(f64::from(header.time)).filter(|time| time.is_finite())
            }
            None => None,
        };
        self.source.file.seek(SeekFrom::Start(offset)).map_err(io)?;
        self.next_offset = offset;
        self.next_frame = target;
        self.source.step = target;
        self.previous_time = previous_time;
        self.last_step = None;
        self.last_precision = None;
        Ok(())
    }

    fn share_index(&mut self) -> Option<crate::FrameIndex> {
        self.frame_index().ok()
    }

    fn adopt_index(&mut self, index: &crate::FrameIndex) -> bool {
        self.adopt_frame_index(index).is_ok()
    }
}
