//! Reading a whole trajectory into memory.

use super::super::{TrajectoryData, TrajectoryFormat, TrajectoryIoError, TrajectoryReadOptions};
use super::{binary, read_bytes};
use binary::{read_dms, read_gsd, read_tng};
use std::path::Path;

/// Reads and fully materialises a trajectory by suffix or format override.
///
/// Memory necessarily scales with the complete output. Use
/// [`super::super::read_trajectory`] for bounded pull-based file processing.
///
/// # Errors
///
/// Returns a dispatch, unit, filesystem or format-specific error without
/// publishing partial frames.
pub fn read_trajectory_materialized(
    path: &Path,
    options: &TrajectoryReadOptions,
) -> Result<TrajectoryData, TrajectoryIoError> {
    let format = options
        .format
        .or_else(|| TrajectoryFormat::infer(path))
        .ok_or(TrajectoryIoError::UnknownFormat)?;
    match format {
        TrajectoryFormat::Xtc => super::super::stream::materialize_xtc(path),
        TrajectoryFormat::Tng => read_tng(path),
        TrajectoryFormat::Gsd => read_gsd(path, options),
        TrajectoryFormat::Dms => read_dms(path),
        _ => read_bytes(path, format, options),
    }
}
