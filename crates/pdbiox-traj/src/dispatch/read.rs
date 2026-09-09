//! Path-based trajectory reading.

use pdbiox_core::io::{InputBuffer, Limits};
use std::path::Path;

use crate::numeric::f32_triplet;

use super::{
    FormatMetadata, TrajectoryData, TrajectoryFormat, TrajectoryIoError, TrajectoryMetadata,
    TrajectoryReadOptions,
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

fn read_amber_netcdf(bytes: &[u8]) -> Result<TrajectoryData, TrajectoryIoError> {
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

fn read_dms(path: &Path) -> Result<TrajectoryData, TrajectoryIoError> {
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

fn read_amber_ascii(
    bytes: &[u8],
    options: &TrajectoryReadOptions,
) -> Result<TrajectoryData, TrajectoryIoError> {
    let text = std::str::from_utf8(bytes).map_err(|_| TrajectoryIoError::InvalidMetadata)?;
    let interpretation = options
        .amber_ascii
        .ok_or(TrajectoryIoError::MissingTopologyOptions)?;
    let title = text
        .lines()
        .next()
        .ok_or(TrajectoryIoError::InvalidMetadata)?
        .into();
    let frames = crate::parse_amber_ascii_trajectory(
        text,
        interpretation.atom_count,
        interpretation.periodic_box,
    )?;
    Ok(TrajectoryData {
        format: TrajectoryFormat::AmberAscii,
        frames,
        metadata: TrajectoryMetadata {
            steps: None,
            format: FormatMetadata::AmberAscii {
                title,
                atom_count: interpretation.atom_count,
                periodic_box: interpretation.periodic_box,
            },
        },
    })
}

fn read_text_output(
    bytes: &[u8],
    format: TrajectoryFormat,
) -> Result<TrajectoryData, TrajectoryIoError> {
    let text = std::str::from_utf8(bytes).map_err(|_| TrajectoryIoError::InvalidMetadata)?;
    let (frames, metadata) = match format {
        TrajectoryFormat::Gamess => {
            let source = crate::parse_gamess_output(text)?;
            let frames = source
                .frames
                .iter()
                .enumerate()
                .map(|(index, frame)| frame.to_timestep(index))
                .collect();
            (frames, FormatMetadata::Gamess(source))
        }
        TrajectoryFormat::LammpsDump => (crate::parse_lammps_dump(text)?, FormatMetadata::None),
        TrajectoryFormat::Gromos11 => {
            let source = crate::parse_gromos11_trc(text)?;
            (source.frames.clone(), FormatMetadata::Gromos11(source))
        }
        _ => return Err(TrajectoryIoError::UnknownFormat),
    };
    Ok(TrajectoryData {
        format,
        frames,
        metadata: TrajectoryMetadata {
            steps: None,
            format: metadata,
        },
    })
}

fn read_charmm(bytes: &[u8]) -> Result<TrajectoryData, TrajectoryIoError> {
    let text = std::str::from_utf8(bytes).map_err(|_| TrajectoryIoError::InvalidMetadata)?;
    let card = crate::parse_charmm_record(text)?;
    Ok(TrajectoryData {
        format: TrajectoryFormat::CharmmCard,
        frames: vec![crate::Timestep {
            positions: card.atoms.iter().map(|atom| atom.position).collect(),
            ..crate::Timestep::default()
        }],
        metadata: TrajectoryMetadata {
            steps: None,
            format: FormatMetadata::CharmmCard(card),
        },
    })
}

fn read_dlpoly(
    bytes: &[u8],
    format: TrajectoryFormat,
) -> Result<TrajectoryData, TrajectoryIoError> {
    let text = std::str::from_utf8(bytes).map_err(|_| TrajectoryIoError::InvalidMetadata)?;
    let (frames, metadata) = match format {
        TrajectoryFormat::DlPolyConfig => {
            let source = crate::parse_dlpoly_config(text)?;
            (
                vec![source.frame.to_timestep(0)],
                FormatMetadata::DlPolyConfig(source),
            )
        }
        TrajectoryFormat::DlPolyHistory => {
            let source = crate::parse_dlpoly_history(text)?;
            let frames = source
                .frames
                .iter()
                .enumerate()
                .map(|(index, frame)| frame.to_timestep(index))
                .collect();
            (frames, FormatMetadata::DlPolyHistory(source))
        }
        _ => return Err(TrajectoryIoError::UnknownFormat),
    };
    Ok(TrajectoryData {
        format,
        frames,
        metadata: TrajectoryMetadata {
            steps: None,
            format: metadata,
        },
    })
}

fn read_txyz(bytes: &[u8]) -> Result<TrajectoryData, TrajectoryIoError> {
    let text = std::str::from_utf8(bytes).map_err(|_| TrajectoryIoError::InvalidMetadata)?;
    let records = crate::parse_txyz_records(text)?;
    let frames = records
        .iter()
        .enumerate()
        .map(|(frame, record)| crate::Timestep {
            frame,
            positions: record.atoms.iter().map(|atom| atom.position).collect(),
            ..crate::Timestep::default()
        })
        .collect();
    Ok(TrajectoryData {
        format: TrajectoryFormat::Txyz,
        frames,
        metadata: TrajectoryMetadata {
            steps: None,
            format: FormatMetadata::Txyz(records),
        },
    })
}

fn read_aims(bytes: &[u8]) -> Result<TrajectoryData, TrajectoryIoError> {
    let text = std::str::from_utf8(bytes).map_err(|_| TrajectoryIoError::InvalidMetadata)?;
    let geometry = crate::parse_aims_geometry(text)?;
    let frame = crate::Timestep {
        positions: geometry.atoms.iter().map(|atom| atom.position).collect(),
        cell: geometry.cell(),
        ..crate::Timestep::default()
    };
    Ok(TrajectoryData {
        format: TrajectoryFormat::Aims,
        frames: vec![frame],
        metadata: TrajectoryMetadata {
            steps: None,
            format: FormatMetadata::Aims(geometry),
        },
    })
}

fn read_xyz(bytes: &[u8]) -> Result<TrajectoryData, TrajectoryIoError> {
    let text = std::str::from_utf8(bytes).map_err(|_| TrajectoryIoError::InvalidMetadata)?;
    let records = crate::parse_xyz(text).ok_or(TrajectoryIoError::InvalidMetadata)?;
    let frames = records
        .iter()
        .enumerate()
        .map(|(frame, record)| crate::Timestep {
            frame,
            positions: record.atoms.iter().map(|atom| atom.position).collect(),
            ..crate::Timestep::default()
        })
        .collect();
    Ok(TrajectoryData {
        format: TrajectoryFormat::Xyz,
        frames,
        metadata: TrajectoryMetadata {
            steps: None,
            format: FormatMetadata::Xyz(records),
        },
    })
}

fn read_gro(bytes: &[u8]) -> Result<TrajectoryData, TrajectoryIoError> {
    let text = std::str::from_utf8(bytes).map_err(|_| TrajectoryIoError::InvalidMetadata)?;
    let records = crate::parse_gro_records(text)?;
    let frames = records
        .iter()
        .enumerate()
        .map(|(frame, record)| crate::Timestep {
            frame,
            positions: record.atoms.iter().map(|atom| atom.position).collect(),
            velocities: record
                .atoms
                .first()
                .and_then(|atom| atom.velocity)
                .map(|_| {
                    record
                        .atoms
                        .iter()
                        .filter_map(|atom| atom.velocity)
                        .collect()
                }),
            ..crate::Timestep::default()
        })
        .collect();
    Ok(TrajectoryData {
        format: TrajectoryFormat::Gro,
        frames,
        metadata: TrajectoryMetadata {
            steps: None,
            format: FormatMetadata::Gro(records),
        },
    })
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

fn read_tng(path: &Path) -> Result<TrajectoryData, TrajectoryIoError> {
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

fn read_gsd(
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
