//! Borrowed-frame contacts, with memory independent of trajectory length.

use crate::{TrajectoryError, TrajectoryReader, run_analysis_stream};
use molframe_core::{AtomSelection, ExecutionContext};
use molframe_spatial::{PairQuery, PeriodicBox, SpatialSearchOptions, count_pairs_within};

/// Emits one unique contact count per frame in source order.
///
/// Coordinates stay in the reusable reader buffer. Each frame's periodic cell
/// is honored when present; absent cells use ordinary Cartesian distances.
/// The spatial index and bounded parallel reductions share the reader account.
///
/// # Errors
///
/// Returns decoder, invalid cutoff or cell, resource, cancellation or sink errors.
pub fn contact_counts_stream<R: TrajectoryReader + ?Sized>(
    reader: &mut R,
    cutoff: f32,
    options: SpatialSearchOptions,
    context: &ExecutionContext,
    frame_workspace_bytes: usize,
    mut emit: impl FnMut(usize, Option<f64>, u64) -> Result<(), TrajectoryError>,
) -> Result<u64, TrajectoryError> {
    let atoms = u32::try_from(reader.n_atoms()).map_err(|_| TrajectoryError::IdentityOverflow)?;
    if !cutoff.is_finite() || cutoff < 0.0 {
        return Err(molframe_spatial::SpatialError::InvalidCutoff.into());
    }
    options.plan(atoms as usize, atoms as usize, cutoff)?;
    let all = AtomSelection::All(atoms);
    run_analysis_stream(reader, context, frame_workspace_bytes, |frame, context| {
        let periodic = frame.cell.map(PeriodicBox::from_cell).transpose()?;
        let count = count_pairs_within(&PairQuery {
            positions: &frame.positions,
            left: &all,
            right: &all,
            cutoff,
            options,
            periodic: periodic.as_ref(),
            context,
        })?;
        emit(frame.frame, frame.time, count)
    })
}

#[cfg(test)]
#[path = "stream_contacts_tests.rs"]
mod tests;
