//! Faithful GROMACS GRO coordinate records in canonical ångström units.

use crate::Frame;
use std::fmt::Write as _;

const NM_TO_ANGSTROM: f32 = 10.0;
const COORD_START: usize = 20;
const COORD_WIDTH: usize = 8;

/// One atom line and all topology fields carried by GRO.
#[derive(Clone, Debug, PartialEq)]
pub struct GroAtom {
    /// Residue number modulo the GRO five-column field.
    pub residue_number: u32,
    /// Residue name.
    pub residue_name: Box<str>,
    /// Atom name.
    pub atom_name: Box<str>,
    /// Atom serial modulo the GRO five-column field.
    pub atom_number: u32,
    /// Cartesian position in ångström.
    pub position: [f32; 3],
    /// Velocity in ångström per picosecond, when declared.
    pub velocity: Option<[f32; 3]>,
}

/// One complete GRO frame.
#[derive(Clone, Debug, PartialEq)]
pub struct GroFrame {
    /// Title line.
    pub title: Box<str>,
    /// Atom records in source order.
    pub atoms: Vec<GroAtom>,
    /// Three or nine GRO box scalars, converted from nm to ångström.
    pub box_values: Vec<f64>,
}

impl GroFrame {
    /// Returns its positions through the common coordinate-only view.
    #[must_use]
    pub fn to_frame(&self) -> Frame {
        Frame {
            positions: self.atoms.iter().map(|atom| atom.position).collect(),
        }
    }
}

/// Invalid GRO syntax or a value outside its fixed-width representation.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum GroError {
    /// A required line or fixed field is malformed.
    #[error("malformed GRO record")]
    Malformed,
    /// Mixed velocity presence or a value overflows a fixed-width field.
    #[error("GRO value cannot be represented faithfully")]
    Unrepresentable,
}

/// Parses complete GRO frames while retaining titles, topology, velocities and boxes.
///
/// # Errors
///
/// Returns [`GroError`] for malformed fixed fields or inconsistent velocity columns.
pub fn parse_gro_records(text: &str) -> Result<Vec<GroFrame>, GroError> {
    let mut frames = Vec::new();
    let mut lines = text.lines();
    loop {
        let Some(title) = lines.next() else {
            return Ok(frames);
        };
        if title.is_empty() && lines.clone().all(|line| line.trim().is_empty()) {
            return Ok(frames);
        }
        let atom_count = number(lines.next())?;
        let mut atoms = Vec::with_capacity(atom_count);
        for _ in 0..atom_count {
            atoms.push(parse_atom(lines.next().ok_or(GroError::Malformed)?)?);
        }
        let box_values = lines
            .next()
            .ok_or(GroError::Malformed)?
            .split_whitespace()
            .map(|value| {
                value
                    .parse::<f64>()
                    .map(|value| value * 10.0)
                    .map_err(|_| GroError::Malformed)
            })
            .collect::<Result<Vec<_>, _>>()?;
        if !matches!(box_values.len(), 3 | 9) {
            return Err(GroError::Malformed);
        }
        frames.push(GroFrame {
            title: title.into(),
            atoms,
            box_values,
        });
    }
}

/// Writes complete GRO records without synthesizing topology or periodic boxes.
///
/// # Errors
///
/// Returns [`GroError::Unrepresentable`] for overflow, invalid box shape, or mixed velocities.
pub fn write_gro(frames: &[GroFrame]) -> Result<String, GroError> {
    let mut output = String::new();
    for frame in frames {
        validate(frame)?;
        writeln!(output, "{}\n{}", frame.title, frame.atoms.len())
            .map_err(|_| GroError::Unrepresentable)?;
        for atom in &frame.atoms {
            write_atom(&mut output, atom)?;
        }
        for value in &frame.box_values {
            write!(output, " {value:>10.5}", value = value / 10.0)
                .map_err(|_| GroError::Unrepresentable)?;
        }
        output.push('\n');
    }
    Ok(output)
}

fn parse_atom(line: &str) -> Result<GroAtom, GroError> {
    if line.len() < 44 {
        return Err(GroError::Malformed);
    }
    let position = [0, 1, 2]
        .map(|axis| fixed_number(line, COORD_START + axis * COORD_WIDTH, 8))
        .into_iter()
        .collect::<Result<Vec<_>, _>>()?;
    let velocity = if line.len() >= 68 {
        Some(
            [0, 1, 2]
                .map(|axis| fixed_number(line, 44 + axis * 8, 8))
                .into_iter()
                .collect::<Result<Vec<_>, _>>()?,
        )
    } else if line.len() == 44 {
        None
    } else {
        return Err(GroError::Malformed);
    };
    Ok(GroAtom {
        residue_number: number(line.get(0..5))?,
        residue_name: line.get(5..10).ok_or(GroError::Malformed)?.trim().into(),
        atom_name: line.get(10..15).ok_or(GroError::Malformed)?.trim().into(),
        atom_number: number(line.get(15..20))?,
        position: [position[0], position[1], position[2]].map(|value| value * NM_TO_ANGSTROM),
        velocity: velocity
            .map(|values| [values[0], values[1], values[2]].map(|value| value * NM_TO_ANGSTROM)),
    })
}

fn validate(frame: &GroFrame) -> Result<(), GroError> {
    let velocity = frame
        .atoms
        .first()
        .is_some_and(|atom| atom.velocity.is_some());
    if frame.title.contains(['\n', '\r'])
        || !matches!(frame.box_values.len(), 3 | 9)
        || frame.box_values.iter().any(|value| !value.is_finite())
        || frame.atoms.iter().any(|atom| {
            atom.residue_number > 99_999
                || atom.atom_number > 99_999
                || atom.residue_name.is_empty()
                || atom.residue_name.len() > 5
                || atom.atom_name.is_empty()
                || atom.atom_name.len() > 5
                || atom.position.iter().any(|value| !value.is_finite())
                || atom.velocity.is_some() != velocity
                || atom
                    .velocity
                    .is_some_and(|values| values.iter().any(|value| !value.is_finite()))
        })
    {
        return Err(GroError::Unrepresentable);
    }
    Ok(())
}

fn write_atom(output: &mut String, atom: &GroAtom) -> Result<(), GroError> {
    write!(
        output,
        "{:>5}{:<5}{:>5}{:>5}",
        atom.residue_number, atom.residue_name, atom.atom_name, atom.atom_number
    )
    .map_err(|_| GroError::Unrepresentable)?;
    for value in atom.position {
        fixed_write(output, value / NM_TO_ANGSTROM, 3)?;
    }
    if let Some(velocity) = atom.velocity {
        for value in velocity {
            fixed_write(output, value / NM_TO_ANGSTROM, 4)?;
        }
    }
    output.push('\n');
    Ok(())
}

fn fixed_write(output: &mut String, value: f32, precision: usize) -> Result<(), GroError> {
    let rendered = format!("{value:>8.precision$}");
    if rendered.len() != 8 {
        return Err(GroError::Unrepresentable);
    }
    output.push_str(&rendered);
    Ok(())
}

fn fixed_number(line: &str, start: usize, width: usize) -> Result<f32, GroError> {
    number(line.get(start..start + width))
}

fn number<T: std::str::FromStr>(value: Option<&str>) -> Result<T, GroError> {
    value
        .and_then(|value| value.trim().parse().ok())
        .ok_or(GroError::Malformed)
}

#[cfg(test)]
#[path = "gro_tests.rs"]
mod tests;
