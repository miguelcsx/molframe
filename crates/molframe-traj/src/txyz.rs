//! Faithful Tinker XYZ/TXYZ and multi-frame ARC records.

use crate::Frame;
use std::collections::BTreeSet;
use std::fmt::Write as _;

/// One Tinker atom with force-field type and connectivity.
#[derive(Clone, Debug, PartialEq)]
pub struct TxyzAtom {
    /// Positive source identifier.
    pub id: u32,
    /// Atom name.
    pub name: Box<str>,
    /// Cartesian position in ångström.
    pub position: [f32; 3],
    /// Force-field atom type retained lexically.
    pub atom_type: Box<str>,
    /// Bonded atom identifiers in source order.
    pub bonds: Vec<u32>,
}

/// One complete Tinker XYZ frame.
#[derive(Clone, Debug, PartialEq)]
pub struct TxyzFrame {
    /// Text following the atom count.
    pub title: Box<str>,
    /// Atoms in source order.
    pub atoms: Vec<TxyzAtom>,
}

impl TxyzFrame {
    /// Returns the coordinate-only common frame.
    #[must_use]
    pub fn to_frame(&self) -> Frame {
        Frame {
            positions: self.atoms.iter().map(|atom| atom.position).collect(),
        }
    }
}

/// Invalid or unrepresentable Tinker XYZ data.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum TxyzError {
    /// Required syntax or connectivity is malformed.
    #[error("malformed Tinker XYZ/ARC record")]
    Malformed,
    /// Text or numeric data cannot be emitted faithfully.
    #[error("Tinker XYZ/ARC value is unrepresentable")]
    Unrepresentable,
}

/// Parses complete TXYZ/ARC records.
///
/// # Errors
///
/// Returns [`TxyzError`] for malformed atoms, duplicate IDs, or dangling bonds.
pub fn parse_txyz_records(text: &str) -> Result<Vec<TxyzFrame>, TxyzError> {
    let mut lines = text.lines();
    let mut frames = Vec::new();
    loop {
        let Some(header) = lines.by_ref().find(|line| !line.trim().is_empty()) else {
            return Ok(frames);
        };
        let mut header_fields = header.trim().splitn(2, char::is_whitespace);
        let count = number(header_fields.next())?;
        let title = header_fields.next().map_or("", str::trim).into();
        let mut atoms = Vec::with_capacity(count);
        let mut ids = BTreeSet::new();
        for _ in 0..count {
            let fields: Vec<_> = lines
                .next()
                .ok_or(TxyzError::Malformed)?
                .split_whitespace()
                .collect();
            if fields.len() < 6 {
                return Err(TxyzError::Malformed);
            }
            let id = number(Some(fields[0]))?;
            if id == 0 || !ids.insert(id) {
                return Err(TxyzError::Malformed);
            }
            atoms.push(TxyzAtom {
                id,
                name: fields[1].into(),
                position: [
                    number(Some(fields[2]))?,
                    number(Some(fields[3]))?,
                    number(Some(fields[4]))?,
                ],
                atom_type: fields[5].into(),
                bonds: fields[6..]
                    .iter()
                    .map(|value| number(Some(value)))
                    .collect::<Result<_, _>>()?,
            });
        }
        if atoms
            .iter()
            .flat_map(|atom| &atom.bonds)
            .any(|bond| !ids.contains(bond))
        {
            return Err(TxyzError::Malformed);
        }
        frames.push(TxyzFrame { title, atoms });
    }
}

/// Writes complete TXYZ/ARC records without generic names, types, or bonds.
///
/// # Errors
///
/// Returns [`TxyzError::Unrepresentable`] for invalid text, IDs, connectivity, or coordinates.
pub fn write_txyz(frames: &[TxyzFrame]) -> Result<String, TxyzError> {
    let mut output = String::new();
    for frame in frames {
        validate(frame)?;
        writeln!(output, "{}  {}", frame.atoms.len(), frame.title)
            .map_err(|_| TxyzError::Unrepresentable)?;
        for atom in &frame.atoms {
            write!(
                output,
                "{} {} {:.7} {:.7} {:.7} {}",
                atom.id,
                atom.name,
                atom.position[0],
                atom.position[1],
                atom.position[2],
                atom.atom_type
            )
            .map_err(|_| TxyzError::Unrepresentable)?;
            for bond in &atom.bonds {
                write!(output, " {bond}").map_err(|_| TxyzError::Unrepresentable)?;
            }
            output.push('\n');
        }
    }
    Ok(output)
}

fn validate(frame: &TxyzFrame) -> Result<(), TxyzError> {
    let ids: BTreeSet<_> = frame.atoms.iter().map(|atom| atom.id).collect();
    if frame.title.contains(['\n', '\r'])
        || ids.len() != frame.atoms.len()
        || frame.atoms.iter().any(|atom| {
            atom.id == 0
                || invalid_token(&atom.name)
                || invalid_token(&atom.atom_type)
                || atom.position.iter().any(|value| !value.is_finite())
                || atom.bonds.iter().any(|bond| !ids.contains(bond))
        })
    {
        Err(TxyzError::Unrepresentable)
    } else {
        Ok(())
    }
}

fn invalid_token(value: &str) -> bool {
    value.is_empty() || value.chars().any(char::is_whitespace)
}

fn number<T: std::str::FromStr>(value: Option<&str>) -> Result<T, TxyzError> {
    value
        .and_then(|value| value.parse().ok())
        .ok_or(TxyzError::Malformed)
}

#[cfg(test)]
#[path = "txyz_tests.rs"]
mod tests;
