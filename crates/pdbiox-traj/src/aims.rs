//! FHI-aims `geometry.in` coordinate reader and writer.

use crate::Frame;
use crate::numeric::f32_triplet;
use pdbiox_core::structure::UnitCell;
use std::fmt::Write;

/// One FHI-aims atom with its species label.
#[derive(Clone, Debug, PartialEq)]
pub struct AimsAtom {
    /// Chemical species as written by FHI-aims.
    pub species: Box<str>,
    /// Cartesian position in ångström.
    pub position: [f32; 3],
}

/// Coordinate-bearing portion of an FHI-aims geometry.
#[derive(Clone, Debug, PartialEq)]
pub struct AimsGeometry {
    /// Atoms in declaration order.
    pub atoms: Vec<AimsAtom>,
    /// Periodic lattice vectors in ångström, when present.
    pub lattice_vectors: Option<[[f64; 3]; 3]>,
}

impl AimsGeometry {
    /// Drops species labels and returns the shared trajectory frame type.
    #[must_use]
    pub fn to_frame(&self) -> Frame {
        Frame {
            positions: self.atoms.iter().map(|atom| atom.position).collect(),
        }
    }

    /// Converts the lattice vectors to crystallographic lengths and angles.
    #[must_use]
    pub fn cell(&self) -> Option<UnitCell> {
        self.lattice_vectors
            .and_then(crate::cell::cell_from_vectors)
    }
}

/// Malformed coordinate-bearing FHI-aims input.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum AimsError {
    /// A recognised coordinate directive has missing fields.
    #[error("incomplete FHI-aims coordinate directive on line {line}")]
    MissingField {
        /// One-based source line.
        line: usize,
    },
    /// One coordinate is not a finite number.
    #[error("invalid FHI-aims coordinate on line {line}")]
    InvalidNumber {
        /// One-based source line.
        line: usize,
    },
    /// A fractional atom occurred without exactly three lattice vectors.
    #[error("fractional FHI-aims atom requires exactly three lattice vectors")]
    FractionalWithoutCell,
    /// More than three lattice vectors were declared.
    #[error("FHI-aims geometry declares more than three lattice vectors")]
    TooManyLatticeVectors,
    /// The lattice vectors do not span a volume.
    #[error("FHI-aims lattice vectors are singular")]
    SingularCell,
}

enum PendingAtom {
    Cartesian(AimsAtom),
    Fractional(Box<str>, [f64; 3], usize),
}

/// Parses Cartesian and fractional atoms plus periodic lattice vectors.
///
/// Non-coordinate FHI-aims directives are intentionally left to the electronic
/// structure input layer and do not alter the coordinate result.
///
/// # Errors
///
/// Returns an error for malformed recognised directives or fractional atoms
/// without a valid three-vector cell.
pub fn parse_aims_geometry(text: &str) -> Result<AimsGeometry, AimsError> {
    let mut vectors = Vec::new();
    let mut atoms = Vec::new();
    for (index, raw) in text.lines().enumerate() {
        let line_number = index + 1;
        let content = match raw.split('#').next() {
            Some(content) => content.trim(),
            None => "",
        };
        if content.is_empty() {
            continue;
        }
        let mut fields = content.split_whitespace();
        match fields.next() {
            Some("lattice_vector") => {
                if vectors.len() == 3 {
                    return Err(AimsError::TooManyLatticeVectors);
                }
                vectors.push(numbers(&mut fields, line_number)?);
            }
            Some("atom") => {
                let position = f32_triplet(numbers(&mut fields, line_number)?)
                    .ok_or(AimsError::InvalidNumber { line: line_number })?;
                let species = species(&mut fields, line_number)?;
                atoms.push(PendingAtom::Cartesian(AimsAtom { species, position }));
            }
            Some("atom_frac") => {
                let position = numbers(&mut fields, line_number)?;
                atoms.push(PendingAtom::Fractional(
                    species(&mut fields, line_number)?,
                    position,
                    line_number,
                ));
            }
            _ => {}
        }
    }
    let lattice_vectors = match vectors.as_slice() {
        [] => None,
        [a, b, c] => Some([*a, *b, *c]),
        _ => return Err(AimsError::FractionalWithoutCell),
    };
    if let Some(vectors) = lattice_vectors
        && crate::cell::determinant(vectors).abs() <= f64::EPSILON
    {
        return Err(AimsError::SingularCell);
    }
    let atoms = atoms
        .into_iter()
        .map(|atom| match atom {
            PendingAtom::Cartesian(atom) => Ok(atom),
            PendingAtom::Fractional(species, position, line) => {
                let vectors = lattice_vectors.ok_or(AimsError::FractionalWithoutCell)?;
                let position = fractional_to_cartesian(vectors, position)
                    .ok_or(AimsError::InvalidNumber { line })?;
                Ok(AimsAtom { species, position })
            }
        })
        .collect::<Result<_, _>>()?;
    Ok(AimsGeometry {
        atoms,
        lattice_vectors,
    })
}

fn numbers<'a>(
    fields: &mut impl Iterator<Item = &'a str>,
    line: usize,
) -> Result<[f64; 3], AimsError> {
    let mut values = [0.0; 3];
    for value in &mut values {
        let raw = fields.next().ok_or(AimsError::MissingField { line })?;
        *value = raw
            .parse()
            .ok()
            .filter(|number: &f64| number.is_finite())
            .ok_or(AimsError::InvalidNumber { line })?;
    }
    Ok(values)
}

fn species<'a>(
    fields: &mut impl Iterator<Item = &'a str>,
    line: usize,
) -> Result<Box<str>, AimsError> {
    fields
        .next()
        .map(Box::from)
        .ok_or(AimsError::MissingField { line })
}

fn fractional_to_cartesian(vectors: [[f64; 3]; 3], value: [f64; 3]) -> Option<[f32; 3]> {
    f32_triplet([0, 1, 2].map(|axis| {
        value
            .iter()
            .zip(vectors.iter())
            .map(|(factor, vector)| factor * vector[axis])
            .sum::<f64>()
    }))
}

/// Writes deterministic Cartesian FHI-aims coordinates and lattice vectors.
#[must_use]
pub fn write_aims_geometry(geometry: &AimsGeometry) -> String {
    let mut output = String::new();
    if let Some(vectors) = geometry.lattice_vectors {
        for vector in vectors {
            let _ = writeln!(
                output,
                "lattice_vector {:.10} {:.10} {:.10}",
                vector[0], vector[1], vector[2]
            );
        }
    }
    for atom in &geometry.atoms {
        let _ = writeln!(
            output,
            "atom {:.10} {:.10} {:.10} {}",
            atom.position[0], atom.position[1], atom.position[2], atom.species
        );
    }
    output
}

#[cfg(test)]
#[path = "aims_tests.rs"]
mod tests;
