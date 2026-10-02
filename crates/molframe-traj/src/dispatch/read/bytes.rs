//! Reading one trajectory from a byte buffer.

use super::super::{
    FormatMetadata, TrajectoryData, TrajectoryFormat, TrajectoryIoError, TrajectoryMetadata,
    TrajectoryReadOptions,
};
use super::{binary, dcd_steps, frame_steps, integer_data, read_snapshot, text};
use binary::read_amber_netcdf;
use molframe_core::io::{InputBuffer, Limits};
use std::path::Path;
use text::{
    read_aims, read_amber_ascii, read_charmm, read_dlpoly, read_gro, read_text_output, read_txyz,
    read_xyz,
};

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
pub(super) fn read_bytes(
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
