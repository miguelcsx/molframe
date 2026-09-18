//! Path-based trajectory reading.
//!
//! This module owns format selection and the shared bounded-input path. The two
//! siblings own the per-format decoders, split by whether the payload is
//! decoded as text or read as raw bytes.

mod binary;
mod text;

use molframe_core::io::{InputBuffer, Limits};
use std::path::Path;

use super::{
    FormatMetadata, TrajectoryData, TrajectoryFormat, TrajectoryIoError, TrajectoryMetadata,
    TrajectoryReadOptions,
};
use binary::{read_amber_netcdf, read_dms, read_gsd, read_tng};
use text::{
    read_aims, read_amber_ascii, read_charmm, read_dlpoly, read_gro, read_text_output, read_txyz,
    read_xyz,
};

/// Reads and fully materialises a trajectory by suffix or format override.
///
/// Memory necessarily scales with the complete output. Use
/// [`super::read_trajectory`] for bounded pull-based file processing.
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
        TrajectoryFormat::Xtc => super::stream::materialize_xtc(path),
        TrajectoryFormat::Tng => read_tng(path),
        TrajectoryFormat::Gsd => read_gsd(path, options),
        TrajectoryFormat::Dms => read_dms(path),
        _ => read_bytes(path, format, options),
    }
}

/// Reads a whole-file format through the shared bounded input abstraction.
///
/// `std::fs::read` would put the file in anonymous memory, which must fit in
/// RAM and swap. Going through [`InputBuffer`] instead gives file-backed pages
/// the kernel can evict once they have been passed, which is the difference
/// between a large trajectory that reads slowly and one that cannot be read at
/// all. It also brings transparent decompression, which the direct read did not
/// have.
///
/// This is not yet streaming: these formats still address the whole input at
/// once, and making them windowed is per-format work. What changes here is the
/// kind of memory the whole input occupies.
fn read_bytes(
    path: &Path,
    format: TrajectoryFormat,
    options: &TrajectoryReadOptions,
) -> Result<TrajectoryData, TrajectoryIoError> {
    let input = InputBuffer::open(path, Limits::default())
        .map_err(|finding| TrajectoryIoError::Input(Box::new(finding)))?;
    let bytes = input.as_bytes();
    let data = match format {
        TrajectoryFormat::Trr => {
            let source = crate::parse_trr(bytes)?;
            TrajectoryData {
                format,
                frames: source.frames,
                metadata: TrajectoryMetadata {
                    steps: Some(source.steps.into_iter().map(i64::from).collect()),
                    format: FormatMetadata::Trr {
                        precision: source.precision,
                    },
                },
            }
        }
        TrajectoryFormat::Dcd => {
            let source = crate::parse_dcd(bytes)?;
            let steps = dcd_steps(&source.header)?;
            TrajectoryData {
                format,
                frames: source.frames,
                metadata: TrajectoryMetadata {
                    steps: Some(steps),
                    format: FormatMetadata::Dcd(source.header),
                },
            }
        }
        TrajectoryFormat::AmberNetcdf => read_amber_netcdf(bytes)?,
        TrajectoryFormat::H5md => {
            let source = crate::parse_h5md_record_with_options(bytes, &options.h5md)?;
            let steps = frame_steps(&source.frames)?;
            TrajectoryData {
                format,
                frames: source.frames,
                metadata: TrajectoryMetadata {
                    steps: Some(steps),
                    format: FormatMetadata::H5md(source.metadata),
                },
            }
        }
        TrajectoryFormat::Trz => {
            let source = crate::parse_trz(bytes)?;
            let steps = source
                .frames
                .iter()
                .map(|frame| integer_data(frame, "trajectory_step"))
                .collect::<Result<Vec<_>, _>>()?;
            TrajectoryData {
                format,
                frames: source.frames,
                metadata: TrajectoryMetadata {
                    steps: Some(steps),
                    format: FormatMetadata::Trz {
                        title: source.title,
                        has_forces: source.has_forces,
                    },
                },
            }
        }
        TrajectoryFormat::Namd | TrajectoryFormat::AmberRestart => {
            read_snapshot(bytes, format, options)?
        }
        TrajectoryFormat::AmberAscii => read_amber_ascii(bytes, options)?,
        TrajectoryFormat::Gro => read_gro(bytes)?,
        TrajectoryFormat::Xyz => read_xyz(bytes)?,
        TrajectoryFormat::Aims => read_aims(bytes)?,
        TrajectoryFormat::Txyz => read_txyz(bytes)?,
        TrajectoryFormat::DlPolyConfig | TrajectoryFormat::DlPolyHistory => {
            read_dlpoly(bytes, format)?
        }
        TrajectoryFormat::CharmmCard => read_charmm(bytes)?,
        TrajectoryFormat::Gamess | TrajectoryFormat::LammpsDump | TrajectoryFormat::Gromos11 => {
            read_text_output(bytes, format)?
        }
        TrajectoryFormat::Xtc
        | TrajectoryFormat::Tng
        | TrajectoryFormat::Gsd
        | TrajectoryFormat::Dms => {
            return Err(TrajectoryIoError::UnknownFormat);
        }
    };
    Ok(data)
}

fn read_snapshot(
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

fn dcd_steps(header: &crate::DcdHeader) -> Result<Vec<i64>, TrajectoryIoError> {
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

fn frame_steps(frames: &[crate::Timestep]) -> Result<Vec<i64>, TrajectoryIoError> {
    frames
        .iter()
        .map(|frame| i64::try_from(frame.frame).map_err(|_| TrajectoryIoError::InvalidMetadata))
        .collect()
}

fn integer_data(frame: &crate::Timestep, key: &str) -> Result<i64, TrajectoryIoError> {
    match frame.data.get(key) {
        Some(crate::FrameValue::Integer(value)) => Ok(*value),
        _ => Err(TrajectoryIoError::InvalidMetadata),
    }
}
