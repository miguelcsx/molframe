//! Reader ownership shared by native, CLI and Python streaming consumers.

use super::{TrajectoryIoError, TrajectoryReaderOptions, read_trajectory};
use crate::{RandomAccess, Timestep, TrajectoryError, TrajectoryReader, Units};
use pdbiox_core::{CancellationToken, ExecutionContext, MemoryReservation};
use std::path::Path;

/// Opens a pull reader after reserving its complete configured allowance.
///
/// The reader retains this charge until dropped. `memory_limit_bytes` includes
/// the ownership wrapper; the decoder receives the remaining allowance. Caller
/// frame buffers and retained outputs require their own reservations, as in
/// [`crate::run_analysis_stream`]. Context clones share the same account.
///
/// # Errors
///
/// Returns a memory error before opening the input if admission fails. Format,
/// capability and decoder errors release the reservation before returning.
pub fn read_trajectory_in(
    path: &Path,
    options: &TrajectoryReaderOptions,
    context: &ExecutionContext,
) -> Result<Box<dyn TrajectoryReader>, TrajectoryIoError> {
    if context.cancellation().is_cancelled() {
        return Err(TrajectoryError::Cancelled.into());
    }
    let bytes = options.memory_limit_bytes;
    let reservation = context
        .try_reserve(bytes)
        .map_err(|_| TrajectoryError::MemoryLimit {
            required: bytes,
            limit: context
                .memory_budget()
                .bytes()
                .saturating_sub(context.reserved_bytes()),
        })?;
    let decoder_bytes = bytes
        .checked_sub(std::mem::size_of::<AccountedReader>())
        .ok_or(TrajectoryError::MemoryLimit {
            required: std::mem::size_of::<AccountedReader>(),
            limit: bytes,
        })?;
    let reader = read_trajectory(
        path,
        &TrajectoryReaderOptions {
            memory_limit_bytes: decoder_bytes,
            ..*options
        },
    )?;
    Ok(Box::new(AccountedReader {
        reader,
        reservation,
        decoder_bytes,
        cancellation: context.cancellation().clone(),
    }))
}

struct AccountedReader {
    reader: Box<dyn TrajectoryReader>,
    reservation: MemoryReservation,
    decoder_bytes: usize,
    cancellation: CancellationToken,
}

impl TrajectoryReader for AccountedReader {
    fn format(&self) -> &'static str {
        self.reader.format()
    }
    fn n_atoms(&self) -> usize {
        self.reader.n_atoms()
    }
    fn n_frames(&self) -> Option<usize> {
        self.reader.n_frames()
    }
    fn units(&self) -> Units {
        self.reader.units()
    }
    fn random_access(&self) -> RandomAccess {
        self.reader.random_access()
    }
    fn read_next(&mut self, frame: &mut Timestep) -> Result<bool, TrajectoryError> {
        self.read_next_bounded(frame, self.decoder_bytes)
    }
    fn read_next_bounded(
        &mut self,
        frame: &mut Timestep,
        bytes: usize,
    ) -> Result<bool, TrajectoryError> {
        if self.cancellation.is_cancelled() {
            return Err(TrajectoryError::Cancelled);
        }
        debug_assert!(self.reservation.bytes() >= self.decoder_bytes);
        self.reader
            .read_next_bounded(frame, bytes.min(self.decoder_bytes))
    }
    fn seek(&mut self, frame: usize) -> Result<(), TrajectoryError> {
        if self.cancellation.is_cancelled() {
            return Err(TrajectoryError::Cancelled);
        }
        self.reader.seek(frame)
    }
}

#[cfg(test)]
#[path = "accounted_tests.rs"]
mod tests;
