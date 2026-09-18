//! Formatted AMBER restart/inpcrd and ASCII mdcrd readers.

use crate::Timestep;
use molframe_core::structure::UnitCell;

/// AMBER internal velocity units to ångström per picosecond.
pub(crate) const AMBER_VELOCITY_TO_ANGSTROM_PER_PS: f32 = 20.455;

/// Explicit interpretation of values following restart coordinates.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum AmberRestartLayout {
    /// Infer from value count and presence of time on the atom-count line.
    #[default]
    Auto,
    /// Coordinates only.
    Coordinates,
    /// Coordinates followed by three box lengths.
    CoordinatesBox3,
    /// Coordinates followed by three lengths and three angles.
    CoordinatesBox6,
    /// Coordinates followed by velocities.
    CoordinatesVelocities,
    /// Coordinates, velocities and three box lengths.
    CoordinatesVelocitiesBox3,
    /// Coordinates, velocities, three lengths and three angles.
    CoordinatesVelocitiesBox6,
}

/// Malformed formatted AMBER content.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum AmberError {
    /// Atom count is absent, zero or malformed.
    #[error("AMBER atom count is absent or invalid")]
    InvalidAtomCount,
    /// A fixed-width numeric field is malformed.
    #[error("invalid AMBER numeric field")]
    InvalidNumber,
    /// Values do not match the requested atom count and layout.
    #[error("AMBER value count mismatch: expected {expected}, found {found}")]
    ValueCountMismatch {
        /// Required scalar count.
        expected: usize,
        /// Parsed scalar count.
        found: usize,
    },
    /// Automatic restart interpretation has two defensible layouts.
    #[error("AMBER restart layout is ambiguous; select it explicitly")]
    AmbiguousLayout,
    /// A canonical field cannot be represented by formatted AMBER restart.
    #[error("value cannot be represented by formatted AMBER restart")]
    Unrepresentable,
}

/// Reads one formatted AMBER restart or inpcrd record.
///
/// Numeric fields are parsed at the specified `12.7` width, so adjacent
/// negative or overflowing whitespace-free fields remain separable. Stored
/// AMBER velocities are converted to canonical ångström per picosecond.
///
/// # Errors
///
/// Returns an explicit count, number or ambiguity error; it never guesses an
/// ambiguous one- or two-atom tail.
pub fn parse_amber_restart(text: &str, layout: AmberRestartLayout) -> Result<Timestep, AmberError> {
    parse_amber_restart_record(text, layout).map(|record| record.timestep)
}

/// Reads one formatted AMBER restart while preserving its title and resolved layout.
///
/// # Errors
///
/// Returns [`AmberError`] for invalid numeric fields, inconsistent counts, or
/// an ambiguous automatic layout.
pub fn parse_amber_restart_record(
    text: &str,
    layout: AmberRestartLayout,
) -> Result<crate::AmberRestart, AmberError> {
    let mut lines = text.lines();
    let title: Box<str> = lines.next().ok_or(AmberError::InvalidAtomCount)?.into();
    let header = lines.next().ok_or(AmberError::InvalidAtomCount)?;
    let header_fields: Vec<_> = header.split_whitespace().collect();
    let atom_count: usize = header_fields
        .first()
        .and_then(|value| value.parse().ok())
        .filter(|count| *count > 0)
        .ok_or(AmberError::InvalidAtomCount)?;
    let time = header_fields.get(1).and_then(|value| value.parse().ok());
    let values = fixed_values(lines, 12)?;
    let resolved = resolve_layout(layout, atom_count, values.len(), time.is_some())?;
    let timestep = restart_timestep(atom_count, time, &values, resolved)?;
    Ok(crate::AmberRestart {
        title,
        layout: resolved,
        timestep,
    })
}

fn resolve_layout(
    layout: AmberRestartLayout,
    atoms: usize,
    found: usize,
    has_time: bool,
) -> Result<AmberRestartLayout, AmberError> {
    if layout != AmberRestartLayout::Auto {
        return Ok(layout);
    }
    let coordinates = atoms * 3;
    if found < coordinates {
        return Err(AmberError::ValueCountMismatch {
            expected: coordinates,
            found,
        });
    }
    let tail = found - coordinates;
    match (has_time, tail) {
        (_, 0) => Ok(AmberRestartLayout::Coordinates),
        (false, 3) if atoms != 1 => Ok(AmberRestartLayout::CoordinatesBox3),
        (false, 6) if atoms != 2 => Ok(AmberRestartLayout::CoordinatesBox6),
        (true, tail) if tail == coordinates => Ok(AmberRestartLayout::CoordinatesVelocities),
        (true, tail) if tail == coordinates + 3 => {
            Ok(AmberRestartLayout::CoordinatesVelocitiesBox3)
        }
        (true, tail) if tail == coordinates + 6 => {
            Ok(AmberRestartLayout::CoordinatesVelocitiesBox6)
        }
        _ => Err(AmberError::AmbiguousLayout),
    }
}

fn restart_timestep(
    atoms: usize,
    time: Option<f64>,
    values: &[f32],
    layout: AmberRestartLayout,
) -> Result<Timestep, AmberError> {
    let (velocities, box_values) = layout_parts(layout);
    let expected = atoms * 3 * (usize::from(velocities) + 1) + box_values;
    if values.len() != expected {
        return Err(AmberError::ValueCountMismatch {
            expected,
            found: values.len(),
        });
    }
    let coordinate_end = atoms * 3;
    let positions = triples(&values[..coordinate_end]);
    let velocity_end = coordinate_end + usize::from(velocities) * atoms * 3;
    let velocity_values = velocities.then(|| {
        triples(&values[coordinate_end..velocity_end])
            .into_iter()
            .map(|velocity| velocity.map(|value| value * AMBER_VELOCITY_TO_ANGSTROM_PER_PS))
            .collect()
    });
    Ok(Timestep {
        time,
        positions,
        velocities: velocity_values,
        cell: cell(&values[velocity_end..]),
        ..Timestep::default()
    })
}

pub(crate) const fn layout_parts(layout: AmberRestartLayout) -> (bool, usize) {
    match layout {
        AmberRestartLayout::Coordinates | AmberRestartLayout::Auto => (false, 0),
        AmberRestartLayout::CoordinatesBox3 => (false, 3),
        AmberRestartLayout::CoordinatesBox6 => (false, 6),
        AmberRestartLayout::CoordinatesVelocities => (true, 0),
        AmberRestartLayout::CoordinatesVelocitiesBox3 => (true, 3),
        AmberRestartLayout::CoordinatesVelocitiesBox6 => (true, 6),
    }
}

/// Reads a formatted AMBER ASCII trajectory with an explicit topology size.
///
/// Classic mdcrd carries no atom count or box flag, so both are caller inputs.
/// Frames use the specified `8.3` width and optionally end in three box lengths.
///
/// # Errors
///
/// Returns a count or numeric error when the scalar stream cannot be divided
/// into complete frames.
pub fn parse_amber_ascii_trajectory(
    text: &str,
    atom_count: usize,
    periodic_box: bool,
) -> Result<Vec<Timestep>, AmberError> {
    if atom_count == 0 {
        return Err(AmberError::InvalidAtomCount);
    }
    let mut lines = text.lines();
    let _title = lines.next().ok_or(AmberError::InvalidNumber)?;
    let values = fixed_values(lines, 8)?;
    let frame_values = atom_count * 3 + if periodic_box { 3 } else { 0 };
    if values.len() % frame_values != 0 {
        return Err(AmberError::ValueCountMismatch {
            expected: frame_values,
            found: values.len() % frame_values,
        });
    }
    Ok(values
        .chunks_exact(frame_values)
        .enumerate()
        .map(|(frame, values)| {
            let coordinate_end = atom_count * 3;
            Timestep {
                frame,
                positions: triples(&values[..coordinate_end]),
                cell: cell(&values[coordinate_end..]),
                ..Timestep::default()
            }
        })
        .collect())
}

fn fixed_values<'a>(
    lines: impl Iterator<Item = &'a str>,
    width: usize,
) -> Result<Vec<f32>, AmberError> {
    let mut values = Vec::new();
    for line in lines {
        let mut start = 0usize;
        while start < line.len() {
            let end = (start + width).min(line.len());
            let field = line
                .get(start..end)
                .ok_or(AmberError::InvalidNumber)?
                .trim();
            if !field.is_empty() {
                values.push(field.parse().map_err(|_| AmberError::InvalidNumber)?);
            }
            start = end;
        }
    }
    Ok(values)
}

fn triples(values: &[f32]) -> Vec<[f32; 3]> {
    let (rows, _) = values.as_chunks::<3>();
    rows.iter()
        .map(|value| [value[0], value[1], value[2]])
        .collect()
}

fn cell(values: &[f32]) -> Option<UnitCell> {
    match values {
        [a, b, c] => Some(UnitCell {
            lengths: [f64::from(*a), f64::from(*b), f64::from(*c)],
            angles: [90.0; 3],
        }),
        [a, b, c, alpha, beta, gamma] => Some(UnitCell {
            lengths: [f64::from(*a), f64::from(*b), f64::from(*c)],
            angles: [f64::from(*alpha), f64::from(*beta), f64::from(*gamma)],
        }),
        _ => None,
    }
}

#[cfg(test)]
#[path = "amber_tests.rs"]
mod tests;
