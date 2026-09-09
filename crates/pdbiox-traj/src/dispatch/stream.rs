//! Bounded pull-reader dispatch and explicit materialisation.

use std::path::Path;

use crate::{DcdReader, Timestep, TrajectoryReader, XtcError, XtcReader};

use super::{
    FormatMetadata, TrajectoryData, TrajectoryFormat, TrajectoryIoError, TrajectoryMetadata,
    TrajectoryReaderOptions,
};

/// Opens a bounded, caller-buffered trajectory source.
///
/// XTC, DCD and TRR are supported without reading the complete file. Other
/// containers return an explicit capability error until they have a
/// format-native pull decoder.
///
/// # Errors
///
/// Returns format inference, capability, filesystem, record or memory errors.
pub fn read_trajectory(
    path: &Path,
    options: &TrajectoryReaderOptions,
) -> Result<Box<dyn TrajectoryReader>, TrajectoryIoError> {
    let format = options
        .format
        .or_else(|| TrajectoryFormat::infer(path))
        .ok_or(TrajectoryIoError::UnknownFormat)?;
    match format {
        TrajectoryFormat::Xtc => Ok(Box::new(XtcReader::open_with_memory_limit(
            path,
            options.memory_limit_bytes,
        )?)),
        TrajectoryFormat::Dcd => Ok(Box::new(DcdReader::open_with_memory_limit(
            path,
            options.memory_limit_bytes,
        )?)),
        TrajectoryFormat::Trr => crate::TrrReader::open(path, options.memory_limit_bytes)
            .map(|reader| Box::new(reader) as Box<dyn TrajectoryReader>)
            .map_err(TrajectoryIoError::from),
        _ => Err(TrajectoryIoError::PullReaderUnavailable { format }),
    }
}

pub(super) fn materialize_xtc(path: &Path) -> Result<TrajectoryData, TrajectoryIoError> {
    let mut reader = XtcReader::open(path)?;
    let mut frame = Timestep::default();
    let mut frames = Vec::new();
    let mut steps = Vec::new();
    let mut precision = Vec::new();
    while reader.read_next_frame(&mut frame)? {
        steps.push(i64::from(reader.last_step().ok_or(XtcError::InvalidFrame)?));
        precision.push(reader.last_precision().ok_or(XtcError::InvalidFrame)?);
        frames.push(std::mem::take(&mut frame));
    }
    if let Some(dt) = frames.get(1).and_then(|frame| frame.dt)
        && let Some(first) = frames.first_mut()
    {
        first.dt = Some(dt);
    }
    Ok(TrajectoryData {
        format: TrajectoryFormat::Xtc,
        frames,
        metadata: TrajectoryMetadata {
            steps: Some(steps),
            format: FormatMetadata::Xtc { precision },
        },
    })
}
