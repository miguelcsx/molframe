//! Deterministic `DL_POLY` text serialisation.

use crate::dlpoly::{DlPolyConfig, DlPolyError, DlPolyFrame, DlPolyHistory, FORCE_TO_CANONICAL};
use std::fmt::Write;

/// Writes a deterministic `DL_POLY` `CONFIG` representation.
///
/// # Errors
///
/// Returns a consistency error when metadata, cell or level-specific arrays do
/// not match the atom count.
pub fn write_dlpoly_config(config: &DlPolyConfig) -> Result<String, DlPolyError> {
    validate_frame(
        &config.frame,
        config.atoms.len(),
        config.level,
        config.boundary,
    )?;
    let mut output = format!(
        "{}\n{:>10}{:>10}{:>10}\n",
        config.title,
        config.level,
        config.boundary,
        config.atoms.len()
    );
    write_cell(&mut output, config.frame.lattice_vectors);
    for (index, atom) in config.atoms.iter().enumerate() {
        let serial = match atom.index {
            Some(serial) => serial,
            None => u32::try_from(index)
                .map_err(|_| DlPolyError::TopologyDrift)?
                .checked_add(1)
                .ok_or(DlPolyError::TopologyDrift)?,
        };
        let _ = writeln!(output, "{:<8}{serial:>10}", atom.name);
        write_vector(&mut output, config.frame.positions[index]);
        if let Some(values) = &config.frame.velocities {
            write_vector(&mut output, values[index]);
        }
        if let Some(values) = &config.frame.forces {
            write_vector(
                &mut output,
                values[index].map(|value| value / FORCE_TO_CANONICAL),
            );
        }
    }
    Ok(output)
}

/// Writes a deterministic formatted `DL_POLY` `HISTORY` trajectory.
///
/// # Errors
///
/// Returns a consistency error before publishing partial output.
pub fn write_dlpoly_history(history: &DlPolyHistory) -> Result<String, DlPolyError> {
    for frame in &history.frames {
        validate_frame(frame, history.atoms.len(), history.level, history.boundary)?;
    }
    let mut output = format!(
        "{}\n{:>10}{:>10}{:>10}\n",
        history.title,
        history.level,
        history.boundary,
        history.atoms.len()
    );
    for frame in &history.frames {
        let _ = writeln!(
            output,
            "timestep{:>10}{:>10}{:>10}{:>10}{:>16.6}",
            frame.step,
            history.atoms.len(),
            history.level,
            history.boundary,
            frame.time
        );
        write_cell(&mut output, frame.lattice_vectors);
        for (index, atom) in history.atoms.iter().enumerate() {
            let serial = match atom.index {
                Some(serial) => serial,
                None => u32::try_from(index)
                    .map_err(|_| DlPolyError::TopologyDrift)?
                    .checked_add(1)
                    .ok_or(DlPolyError::TopologyDrift)?,
            };
            let mass = atom.mass.ok_or(DlPolyError::TopologyDrift)?;
            let charge = atom.charge.ok_or(DlPolyError::TopologyDrift)?;
            let _ = writeln!(
                output,
                "{:<8}{serial:>10}{mass:>16.8}{charge:>16.8}",
                atom.name
            );
            write_vector(&mut output, frame.positions[index]);
            if let Some(values) = &frame.velocities {
                write_vector(&mut output, values[index]);
            }
            if let Some(values) = &frame.forces {
                write_vector(
                    &mut output,
                    values[index].map(|value| value / FORCE_TO_CANONICAL),
                );
            }
        }
    }
    Ok(output)
}

fn validate_frame(
    frame: &DlPolyFrame,
    atoms: usize,
    level: u8,
    boundary: i32,
) -> Result<(), DlPolyError> {
    let counts_match = frame.positions.len() == atoms
        && frame.velocities.as_ref().map(Vec::len) == (level >= 1).then_some(atoms)
        && frame.forces.as_ref().map(Vec::len) == (level >= 2).then_some(atoms);
    if !counts_match {
        return Err(DlPolyError::CountMismatch {
            expected: atoms,
            found: frame.positions.len(),
        });
    }
    if (boundary == 0) != frame.lattice_vectors.is_none()
        || frame
            .lattice_vectors
            .is_some_and(|vectors| crate::cell::determinant(vectors).abs() <= f64::EPSILON)
    {
        return Err(DlPolyError::InvalidCell);
    }
    Ok(())
}

fn write_cell(output: &mut String, vectors: Option<[[f64; 3]; 3]>) {
    if let Some(vectors) = vectors {
        for vector in vectors {
            let _ = writeln!(
                output,
                "{:>20.10}{:>20.10}{:>20.10}",
                vector[0], vector[1], vector[2]
            );
        }
    }
}

fn write_vector(output: &mut String, vector: [f32; 3]) {
    let _ = writeln!(
        output,
        "{:>20.10}{:>20.10}{:>20.10}",
        vector[0], vector[1], vector[2]
    );
}
