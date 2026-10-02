//! Reading a single-frame snapshot file.

use super::super::{
    FormatMetadata, TrajectoryData, TrajectoryFormat, TrajectoryIoError, TrajectoryMetadata,
    TrajectoryReadOptions,
};

pub(super) fn read_snapshot(
    bytes: &[u8],
    format: TrajectoryFormat,
    options: &TrajectoryReadOptions,
) -> Result<TrajectoryData, TrajectoryIoError> {
    let (frame, metadata) = match format {
        TrajectoryFormat::Namd => {
            let source = crate::parse_namd_binary(bytes)?;
            (source.frame, FormatMetadata::Namd(source.endian))
        }
        TrajectoryFormat::AmberRestart => {
            let text =
                std::str::from_utf8(bytes).map_err(|_| TrajectoryIoError::InvalidMetadata)?;
            let source = crate::parse_amber_restart_record(text, options.amber_restart_layout)?;
            (
                source.timestep,
                FormatMetadata::AmberRestart {
                    title: source.title,
                    layout: source.layout,
                },
            )
        }
        _ => return Err(TrajectoryIoError::UnknownFormat),
    };
    Ok(TrajectoryData {
        format,
        frames: vec![frame],
        metadata: TrajectoryMetadata {
            steps: None,
            format: metadata,
        },
    })
}
