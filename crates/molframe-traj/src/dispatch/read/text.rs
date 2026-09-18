//! Decoders for formats whose payload is UTF-8 text.
//!
//! Every entry point here decodes the whole input as one `&str` and hands it to
//! the format's parser. A payload that is not valid UTF-8 is a metadata error,
//! not a partial read: these files are line-oriented, so a bad byte early on
//! makes every later line unreachable.

use super::super::{
    FormatMetadata, TrajectoryData, TrajectoryFormat, TrajectoryIoError, TrajectoryMetadata,
    TrajectoryReadOptions,
};

pub(super) fn read_amber_ascii(
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

pub(super) fn read_text_output(
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

pub(super) fn read_charmm(bytes: &[u8]) -> Result<TrajectoryData, TrajectoryIoError> {
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

pub(super) fn read_dlpoly(
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

pub(super) fn read_txyz(bytes: &[u8]) -> Result<TrajectoryData, TrajectoryIoError> {
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

pub(super) fn read_aims(bytes: &[u8]) -> Result<TrajectoryData, TrajectoryIoError> {
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

pub(super) fn read_xyz(bytes: &[u8]) -> Result<TrajectoryData, TrajectoryIoError> {
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

pub(super) fn read_gro(bytes: &[u8]) -> Result<TrajectoryData, TrajectoryIoError> {
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
