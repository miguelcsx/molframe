//! Classifying fixed-column models before choosing their coordinate storage.

use super::lines::Lines;
use crate::fixed;
use pdbiox_core::diagnostic::{Code, Diagnostic};

/// A model selected by its deposition order and number.
#[derive(Clone, Copy)]
pub(super) struct ModelSpec {
    pub(super) ordinal: usize,
    pub(super) number: i32,
}

struct Candidate<'a> {
    spec: ModelSpec,
    atoms: Vec<&'a str>,
}

/// Models that require independent structures, or `None` for dense storage.
pub(super) fn ragged_models(
    text: &str,
    only_first: bool,
) -> Result<Option<Vec<ModelSpec>>, Diagnostic> {
    if only_first {
        return Ok(None);
    }
    let candidates = candidates(text)?;
    let Some(first) = candidates.first() else {
        return Ok(None);
    };
    if candidates.len() < 2
        || candidates
            .iter()
            .skip(1)
            .all(|candidate| same_atoms(first, candidate))
    {
        return Ok(None);
    }
    Ok(Some(
        candidates.iter().map(|candidate| candidate.spec).collect(),
    ))
}

fn candidates(text: &str) -> Result<Vec<Candidate<'_>>, Diagnostic> {
    let mut candidates = Vec::new();
    let mut current = None;
    for line in Lines::new(text) {
        let line = line?;
        match fixed::record(line.text) {
            "MODEL" => {
                let ordinal = candidates.len();
                let number = fixed::integer(line.text, 11, 14)
                    .and_then(|value| i32::try_from(value).ok())
                    .ok_or_else(|| {
                        Diagnostic::new(Code::E1202)
                            .with_message("MODEL serial is absent or outside the supported range")
                    })?;
                candidates.push(Candidate {
                    spec: ModelSpec { ordinal, number },
                    atoms: Vec::new(),
                });
                current = Some(ordinal);
            }
            "ATOM" | "HETATM" => {
                if current.is_none() {
                    candidates.push(Candidate {
                        spec: ModelSpec {
                            ordinal: 0,
                            number: 1,
                        },
                        atoms: Vec::new(),
                    });
                    current = Some(0);
                }
                if let Some(candidate) = current.and_then(|index| candidates.get_mut(index)) {
                    candidate.atoms.push(line.text);
                }
            }
            "ENDMDL" => current = None,
            _ => {}
        }
    }
    Ok(candidates)
}

fn same_atoms(left: &Candidate<'_>, right: &Candidate<'_>) -> bool {
    left.atoms.len() == right.atoms.len()
        && left
            .atoms
            .iter()
            .zip(&right.atoms)
            .all(|(left, right)| atom_annotation(left) == atom_annotation(right))
}

fn atom_annotation(line: &str) -> (&str, &str) {
    (fixed::raw(line, 1, 30), fixed::raw(line, 55, line.len()))
}

#[cfg(test)]
#[path = "ensemble_tests.rs"]
mod tests;
