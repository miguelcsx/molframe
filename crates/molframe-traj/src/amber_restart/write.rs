//! Fixed-width AMBER restart writing with explicit canonical-unit conversion.

use super::AmberRestart;
use crate::amber::{
    AMBER_VELOCITY_TO_ANGSTROM_PER_PS, AmberError, AmberRestartLayout, layout_parts,
};
use std::fmt::Write as _;

/// Writes a formatted AMBER restart/inpcrd record.
///
/// Positions and cell lengths are supplied in ångström, time in picoseconds,
/// and velocities in ångström per picosecond. Velocities are converted to the
/// AMBER restart unit before serialization.
///
/// # Errors
///
/// Returns [`AmberError::Unrepresentable`] when the selected layout disagrees
/// with the frame or a formatted field would overflow its fixed width.
pub fn write_amber_restart(record: &AmberRestart) -> Result<String, AmberError> {
    validate(record)?;
    let mut output = String::new();
    writeln!(output, "{}", record.title).map_err(|_| AmberError::Unrepresentable)?;
    write!(output, "{:>6}", record.timestep.positions.len())
        .map_err(|_| AmberError::Unrepresentable)?;
    if let Some(time) = record.timestep.time {
        let rendered = format!("{time:>15.7}");
        if rendered.len() > 15 {
            return Err(AmberError::Unrepresentable);
        }
        output.push_str(&rendered);
    }
    output.push('\n');
    let mut values = Vec::with_capacity(expected_scalar_count(record));
    extend_triples(&mut values, &record.timestep.positions);
    if let Some(velocities) = &record.timestep.velocities {
        values.extend(
            velocities
                .iter()
                .flatten()
                .map(|value| f64::from(*value / AMBER_VELOCITY_TO_ANGSTROM_PER_PS)),
        );
    }
    if let Some(cell) = record.timestep.cell {
        values.extend(cell.lengths);
        if matches!(
            record.layout,
            AmberRestartLayout::CoordinatesBox6 | AmberRestartLayout::CoordinatesVelocitiesBox6
        ) {
            values.extend(cell.angles);
        }
    }
    write_values(&values, &mut output)?;
    Ok(output)
}

fn validate(record: &AmberRestart) -> Result<(), AmberError> {
    if record.layout == AmberRestartLayout::Auto
        || record.title.contains(['\n', '\r'])
        || record.timestep.positions.is_empty()
        || record.timestep.positions.len() > 999_999
        || record.timestep.frame != 0
        || record.timestep.dt.is_some()
        || record.timestep.forces.is_some()
        || !record.timestep.data.is_empty()
        || record.timestep.time.is_some_and(|value| !value.is_finite())
        || record
            .timestep
            .positions
            .iter()
            .flatten()
            .any(|value| !value.is_finite())
    {
        return Err(AmberError::Unrepresentable);
    }
    let (needs_velocities, box_values) = layout_parts(record.layout);
    match (&record.timestep.velocities, needs_velocities) {
        (Some(values), true)
            if values.len() == record.timestep.positions.len()
                && values.iter().flatten().all(|value| value.is_finite()) => {}
        (None, false) => {}
        _ => return Err(AmberError::Unrepresentable),
    }
    validate_cell(record, box_values)
}

fn validate_cell(record: &AmberRestart, box_values: usize) -> Result<(), AmberError> {
    match (record.timestep.cell, box_values) {
        (None, 0) => Ok(()),
        (Some(cell), 3 | 6)
            if cell
                .lengths
                .iter()
                .chain(&cell.angles)
                .all(|value| value.is_finite())
                && (box_values == 6
                    || cell
                        .angles
                        .iter()
                        .all(|angle| (*angle - 90.0).abs() < f64::EPSILON)) =>
        {
            Ok(())
        }
        _ => Err(AmberError::Unrepresentable),
    }
}

fn expected_scalar_count(record: &AmberRestart) -> usize {
    let (velocities, box_values) = layout_parts(record.layout);
    record.timestep.positions.len() * 3 * (usize::from(velocities) + 1) + box_values
}

fn extend_triples(output: &mut Vec<f64>, values: &[[f32; 3]]) {
    output.extend(values.iter().flatten().copied().map(f64::from));
}

fn write_values(values: &[f64], output: &mut String) -> Result<(), AmberError> {
    for (index, value) in values.iter().enumerate() {
        let rendered = format!("{value:>12.7}");
        if rendered.len() > 12 {
            return Err(AmberError::Unrepresentable);
        }
        output.push_str(&rendered);
        if (index + 1) % 6 == 0 {
            output.push('\n');
        }
    }
    if !output.ends_with('\n') {
        output.push('\n');
    }
    Ok(())
}
