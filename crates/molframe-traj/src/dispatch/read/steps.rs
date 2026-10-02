//! Step numbers recovered from headers and frames.

use super::super::TrajectoryIoError;

pub(super) fn dcd_steps(header: &crate::DcdHeader) -> Result<Vec<i64>, TrajectoryIoError> {
    (0..header.frame_count)
        .map(|index| {
            i64::from(header.start_step)
                .checked_add(
                    i64::try_from(index).map_err(|_| TrajectoryIoError::InvalidMetadata)?
                        * i64::from(header.save_interval),
                )
                .ok_or(TrajectoryIoError::InvalidMetadata)
        })
        .collect()
}

pub(super) fn frame_steps(frames: &[crate::Timestep]) -> Result<Vec<i64>, TrajectoryIoError> {
    frames
        .iter()
        .map(|frame| i64::try_from(frame.frame).map_err(|_| TrajectoryIoError::InvalidMetadata))
        .collect()
}

pub(super) fn integer_data(frame: &crate::Timestep, key: &str) -> Result<i64, TrajectoryIoError> {
    match frame.data.get(key) {
        Some(crate::FrameValue::Integer(value)) => Ok(*value),
        _ => Err(TrajectoryIoError::InvalidMetadata),
    }
}
