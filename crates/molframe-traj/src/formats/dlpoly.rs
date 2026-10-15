//! `DL_POLY` `CONFIG` and multi-frame `HISTORY` text formats.

use crate::{FrameValue, Timestep};
use std::collections::BTreeMap;

pub(crate) const FORCE_TO_CANONICAL: f32 = 0.01;

/// One `DL_POLY` particle identity and physical metadata.
#[derive(Clone, Debug, PartialEq)]
pub struct DlPolyAtom {
    /// Case-sensitive particle name.
    pub name: Box<str>,
    /// One-based file index, when declared.
    pub index: Option<u32>,
    /// Mass in daltons, present in HISTORY.
    pub mass: Option<f64>,
    /// Charge in proton-charge units, present in HISTORY.
    pub charge: Option<f64>,
}

/// One coordinate state with canonical molframe units.
#[derive(Clone, Debug, PartialEq)]
pub struct DlPolyFrame {
    /// Simulation step.
    pub step: i64,
    /// Simulation time in picoseconds.
    pub time: f64,
    /// Positions in ångström.
    pub positions: Vec<[f32; 3]>,
    /// Velocities in ångström per picosecond.
    pub velocities: Option<Vec<[f32; 3]>>,
    /// Forces in kJ mol⁻¹ Å⁻¹.
    pub forces: Option<Vec<[f32; 3]>>,
    /// Cell vectors as rows, in ångström.
    pub lattice_vectors: Option<[[f64; 3]; 3]>,
}

impl DlPolyFrame {
    /// Converts to the common reusable timestep representation.
    #[must_use]
    pub fn to_timestep(&self, frame: usize) -> Timestep {
        let mut data = BTreeMap::new();
        data.insert("dlpoly.step".into(), FrameValue::Integer(self.step));
        Timestep {
            frame,
            time: Some(self.time),
            positions: self.positions.clone(),
            velocities: self.velocities.clone(),
            forces: self.forces.clone(),
            cell: self
                .lattice_vectors
                .and_then(crate::cell::cell_from_vectors),
            data,
            ..Timestep::default()
        }
    }
}

/// One `DL_POLY` `CONFIG` state.
#[derive(Clone, Debug, PartialEq)]
pub struct DlPolyConfig {
    /// Title record.
    pub title: Box<str>,
    /// Data level: 0 positions, 1 plus velocities, 2 plus forces.
    pub level: u8,
    /// `DL_POLY` periodic boundary key.
    pub boundary: i32,
    /// Particle identities in file order.
    pub atoms: Vec<DlPolyAtom>,
    /// The single coordinate state.
    pub frame: DlPolyFrame,
}

/// A `DL_POLY` `HISTORY` trajectory over fixed particle identities.
#[derive(Clone, Debug, PartialEq)]
pub struct DlPolyHistory {
    /// Simulation title.
    pub title: Box<str>,
    /// Data level shared by every frame.
    pub level: u8,
    /// Periodic boundary key shared by every frame.
    pub boundary: i32,
    /// Stable particle topology.
    pub atoms: Vec<DlPolyAtom>,
    /// Coordinate frames in source order.
    pub frames: Vec<DlPolyFrame>,
}

/// Malformed or inconsistent `DL_POLY` text data.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum DlPolyError {
    /// Title or control record is absent.
    #[error("DL_POLY header is absent or incomplete")]
    MissingHeader,
    /// Data level is outside 0..=2.
    #[error("DL_POLY data level must be 0, 1 or 2")]
    InvalidLevel,
    /// A recognised numeric record is malformed.
    #[error("invalid DL_POLY numeric record on line {line}")]
    InvalidNumber {
        /// One-based source line.
        line: usize,
    },
    /// The declared particle or value count was not met.
    #[error("DL_POLY expected {expected} records, found {found}")]
    CountMismatch {
        /// Required count.
        expected: usize,
        /// Observed count.
        found: usize,
    },
    /// Atom identity or frame controls changed within HISTORY.
    #[error("DL_POLY HISTORY topology or controls changed between frames")]
    TopologyDrift,
    /// A periodic frame has missing or singular cell vectors.
    #[error("DL_POLY periodic cell is missing or singular")]
    InvalidCell,
}

/// Parses a `DL_POLY` `CONFIG` file, including optional velocities and forces.
///
/// # Errors
///
/// Returns an explicit header, count, number, level or cell error.
pub fn parse_dlpoly_config(text: &str) -> Result<DlPolyConfig, DlPolyError> {
    let lines = nonempty_lines(text);
    let title = lines.first().ok_or(DlPolyError::MissingHeader)?.1.into();
    let (header_line, header) = lines.get(1).ok_or(DlPolyError::MissingHeader)?;
    let fields: Vec<_> = header.split_whitespace().collect();
    let level = parse_level(fields.first().copied(), *header_line)?;
    let boundary: i32 = parse(fields.get(1).copied(), *header_line)?;
    let declared = fields.get(2).and_then(|value| value.parse().ok());
    let mut cursor = 2;
    let lattice_vectors = read_cell(&lines, &mut cursor, boundary)?;
    let stride = 2 + usize::from(level);
    let remaining = lines.len() - cursor;
    if !remaining.is_multiple_of(stride) {
        return Err(DlPolyError::CountMismatch {
            expected: remaining.div_ceil(stride) * stride,
            found: remaining,
        });
    }
    let atom_count = match declared {
        Some(count) => count,
        None => remaining / stride,
    };
    if remaining / stride != atom_count {
        return Err(DlPolyError::CountMismatch {
            expected: atom_count,
            found: remaining / stride,
        });
    }
    let mut atoms = Vec::with_capacity(atom_count);
    let mut positions = Vec::with_capacity(atom_count);
    let mut velocities = (level >= 1).then(|| Vec::with_capacity(atom_count));
    let mut forces = (level >= 2).then(|| Vec::with_capacity(atom_count));
    for _ in 0..atom_count {
        let (line_number, atom_line) = lines[cursor];
        cursor += 1;
        let mut atom_fields = atom_line.split_whitespace();
        let name = atom_fields.next().ok_or(DlPolyError::MissingHeader)?.into();
        let index = atom_fields.next().and_then(|value| value.parse().ok());
        atoms.push(DlPolyAtom {
            name,
            index,
            mass: None,
            charge: None,
        });
        positions.push(read_vector(&lines, &mut cursor)?);
        if let Some(values) = &mut velocities {
            values.push(read_vector(&lines, &mut cursor)?);
        }
        if let Some(values) = &mut forces {
            values.push(read_vector(&lines, &mut cursor)?.map(|value| value * FORCE_TO_CANONICAL));
        }
        if atom_line.trim().is_empty() {
            return Err(DlPolyError::InvalidNumber { line: line_number });
        }
    }
    Ok(DlPolyConfig {
        title,
        level,
        boundary,
        atoms,
        frame: DlPolyFrame {
            step: 0,
            time: 0.0,
            positions,
            velocities,
            forces,
            lattice_vectors,
        },
    })
}

/// Parses every frame in a formatted `DL_POLY` `HISTORY` file.
///
/// # Errors
///
/// Rejects truncated frames and any atom-count, control or identity drift.
pub fn parse_dlpoly_history(text: &str) -> Result<DlPolyHistory, DlPolyError> {
    let lines = nonempty_lines(text);
    let title = lines.first().ok_or(DlPolyError::MissingHeader)?.1.into();
    let (header_line, header) = lines.get(1).ok_or(DlPolyError::MissingHeader)?;
    let fields: Vec<_> = header.split_whitespace().collect();
    let level = parse_level(fields.first().copied(), *header_line)?;
    let boundary: i32 = parse(fields.get(1).copied(), *header_line)?;
    let atom_count: usize = parse(fields.get(2).copied(), *header_line)?;
    let mut cursor = 2;
    let mut topology: Option<Vec<DlPolyAtom>> = None;
    let mut frames = Vec::new();
    while cursor < lines.len() {
        let (frame_line, frame_header) = lines[cursor];
        cursor += 1;
        let fields: Vec<_> = frame_header.split_whitespace().collect();
        if fields.first().copied() != Some("timestep") || fields.len() < 6 {
            return Err(DlPolyError::MissingHeader);
        }
        let step = parse(fields.get(1).copied(), frame_line)?;
        let frame_atoms: usize = parse(fields.get(2).copied(), frame_line)?;
        let frame_level = parse_level(fields.get(3).copied(), frame_line)?;
        let frame_boundary: i32 = parse(fields.get(4).copied(), frame_line)?;
        let time = parse(fields.get(5).copied(), frame_line)?;
        if frame_atoms != atom_count || frame_level != level || frame_boundary != boundary {
            return Err(DlPolyError::TopologyDrift);
        }
        let lattice_vectors = read_cell(&lines, &mut cursor, boundary)?;
        let (atoms, positions, velocities, forces) =
            read_history_atoms(&lines, &mut cursor, atom_count, level)?;
        if topology.as_ref().is_some_and(|expected| expected != &atoms) {
            return Err(DlPolyError::TopologyDrift);
        }
        if topology.is_none() {
            topology = Some(atoms);
        }
        frames.push(DlPolyFrame {
            step,
            time,
            positions,
            velocities,
            forces,
            lattice_vectors,
        });
    }
    Ok(DlPolyHistory {
        title,
        level,
        boundary,
        atoms: topology_or_empty(topology),
        frames,
    })
}

fn topology_or_empty(topology: Option<Vec<DlPolyAtom>>) -> Vec<DlPolyAtom> {
    let Some(atoms) = topology else {
        return Vec::new();
    };
    atoms
}

type AtomValues = (
    Vec<DlPolyAtom>,
    Vec<[f32; 3]>,
    Option<Vec<[f32; 3]>>,
    Option<Vec<[f32; 3]>>,
);

fn read_history_atoms(
    lines: &[(usize, &str)],
    cursor: &mut usize,
    atom_count: usize,
    level: u8,
) -> Result<AtomValues, DlPolyError> {
    let mut atoms = Vec::with_capacity(atom_count);
    let mut positions = Vec::with_capacity(atom_count);
    let mut velocities = (level >= 1).then(|| Vec::with_capacity(atom_count));
    let mut forces = (level >= 2).then(|| Vec::with_capacity(atom_count));
    for _ in 0..atom_count {
        let (line, header) = *lines.get(*cursor).ok_or(DlPolyError::CountMismatch {
            expected: atom_count,
            found: atoms.len(),
        })?;
        *cursor += 1;
        let fields: Vec<_> = header.split_whitespace().collect();
        atoms.push(DlPolyAtom {
            name: fields
                .first()
                .ok_or(DlPolyError::MissingHeader)?
                .to_string()
                .into(),
            index: Some(parse(fields.get(1).copied(), line)?),
            mass: Some(parse(fields.get(2).copied(), line)?),
            charge: Some(parse(fields.get(3).copied(), line)?),
        });
        positions.push(read_vector(lines, cursor)?);
        if let Some(values) = &mut velocities {
            values.push(read_vector(lines, cursor)?);
        }
        if let Some(values) = &mut forces {
            values.push(read_vector(lines, cursor)?.map(|value| value * FORCE_TO_CANONICAL));
        }
    }
    Ok((atoms, positions, velocities, forces))
}

fn nonempty_lines(text: &str) -> Vec<(usize, &str)> {
    text.lines()
        .enumerate()
        .filter(|(_, line)| !line.trim().is_empty())
        .map(|(index, line)| (index + 1, line))
        .collect()
}

fn parse_level(value: Option<&str>, line: usize) -> Result<u8, DlPolyError> {
    let level: u8 = parse(value, line)?;
    (level <= 2)
        .then_some(level)
        .ok_or(DlPolyError::InvalidLevel)
}

fn parse<T: std::str::FromStr>(value: Option<&str>, line: usize) -> Result<T, DlPolyError> {
    value
        .and_then(|value| value.parse().ok())
        .ok_or(DlPolyError::InvalidNumber { line })
}

fn read_vector(lines: &[(usize, &str)], cursor: &mut usize) -> Result<[f32; 3], DlPolyError> {
    let (line, record) = *lines.get(*cursor).ok_or(DlPolyError::CountMismatch {
        expected: *cursor + 1,
        found: lines.len(),
    })?;
    *cursor += 1;
    let mut fields = record.split_whitespace();
    Ok([
        parse(fields.next(), line)?,
        parse(fields.next(), line)?,
        parse(fields.next(), line)?,
    ])
}

fn read_cell(
    lines: &[(usize, &str)],
    cursor: &mut usize,
    boundary: i32,
) -> Result<Option<[[f64; 3]; 3]>, DlPolyError> {
    if boundary == 0 {
        return Ok(None);
    }
    let vectors = [
        read_vector(lines, cursor)?.map(f64::from),
        read_vector(lines, cursor)?.map(f64::from),
        read_vector(lines, cursor)?.map(f64::from),
    ];
    if crate::cell::determinant(vectors).abs() <= f64::EPSILON {
        return Err(DlPolyError::InvalidCell);
    }
    Ok(Some(vectors))
}

#[cfg(test)]
#[path = "dlpoly_tests.rs"]
mod tests;
