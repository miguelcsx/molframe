//! Decoders for formats read as raw bytes or through their own file reader.
//!
//! These payloads carry their structure in a binary header, so they are never
//! decoded as text and two of them — TNG and GSD — must seek rather than take a
//! whole in-memory buffer.

use std::path::Path;

use super::super::{
    FormatMetadata, TrajectoryData, TrajectoryFormat, TrajectoryIoError, TrajectoryMetadata,
    TrajectoryReadOptions,
};
use crate::numeric::f32_triplet;

pub(super) fn read_amber_netcdf(bytes: &[u8]) -> Result<TrajectoryData, TrajectoryIoError> {
    let source = crate::parse_amber_netcdf_record(bytes)?;
    Ok(TrajectoryData {
        format: TrajectoryFormat::AmberNetcdf,
        frames: source.frames,
        metadata: TrajectoryMetadata {
            steps: None,
            format: FormatMetadata::AmberNetcdf(source.metadata),
        },
    })
}

pub(super) fn read_dms(path: &Path) -> Result<TrajectoryData, TrajectoryIoError> {
    let source = crate::read_dms(path)?;
    let velocities = source
        .frame
        .velocities
        .iter()
        .copied()
        .collect::<Option<Vec<_>>>()
        .map(|values| {
            values
                .into_iter()
                .map(|value| f32_triplet(value).ok_or(TrajectoryIoError::UnrepresentableValue))
                .collect::<Result<Vec<_>, _>>()
        })
        .transpose()?;
    let frame = crate::Timestep {
        positions: source
            .frame
            .positions
            .iter()
            .map(|position| f32_triplet(*position).ok_or(TrajectoryIoError::UnrepresentableValue))
            .collect::<Result<Vec<_>, _>>()?,
        velocities,
        cell: source
            .frame
            .cell
            .and_then(|cell| crate::cell::cell_from_vectors(cell.vectors)),
        ..crate::Timestep::default()
    };
    Ok(TrajectoryData {
        format: TrajectoryFormat::Dms,
        frames: vec![frame],
        metadata: TrajectoryMetadata {
            steps: None,
            format: FormatMetadata::Dms(Box::new(source)),
        },
    })
}

pub(super) fn read_tng(path: &Path) -> Result<TrajectoryData, TrajectoryIoError> {
    let source = crate::parse_tng(path)?;
    Ok(TrajectoryData {
        format: TrajectoryFormat::Tng,
        frames: source.frames,
        metadata: TrajectoryMetadata {
            steps: Some(source.steps),
            format: FormatMetadata::Tng {
                distance_unit_exponent: source.distance_unit_exponent,
                compression_precision: source.compression_precision,
                compression: source.compression,
            },
        },
    })
}

pub(super) fn read_gsd(
    path: &Path,
    options: &TrajectoryReadOptions,
) -> Result<TrajectoryData, TrajectoryIoError> {
    let units = options.gsd.ok_or(TrajectoryIoError::MissingUnitOptions)?;
    let source = crate::parse_gsd(path, units)?;
    let steps = source
        .steps
        .into_iter()
        .map(|step| i64::try_from(step).map_err(|_| TrajectoryIoError::InvalidMetadata))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(TrajectoryData {
        format: TrajectoryFormat::Gsd,
        frames: source.frames,
        metadata: TrajectoryMetadata {
            steps: Some(steps),
            format: FormatMetadata::Gsd(units),
        },
    })
}
