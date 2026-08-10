//! Occupancy and B-factor sanity checks.
//!
//! These columns carry values that only make sense in a range: an occupancy is a
//! fraction in (0, 1], and a B-factor is a non-negative displacement. A value
//! outside its range is a data problem — a truncated field, a placeholder left
//! in, an occupancy that should have been split across altlocs — and each is
//! reported against the atom it belongs to rather than summarised away.

use pdbiox_core::index::AtomIndex;
use pdbiox_core::structure::Structure;

/// What is wrong with an atom's occupancy or B-factor.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum QualityIssue {
    /// Occupancy is exactly zero, so the atom contributes nothing yet is present.
    ZeroOccupancy,
    /// Occupancy is negative or above one, which no fraction can be.
    OccupancyOutOfRange,
    /// The B-factor is negative, which a mean-square displacement cannot be.
    NegativeBFactor,
}

/// An atom flagged by a quality check, and why.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QualityFlag {
    /// The atom at fault.
    pub atom: AtomIndex,
    /// The problem found.
    pub issue: QualityIssue,
}

/// Flags atoms with an out-of-range occupancy or a negative B-factor.
///
/// An atom can raise more than one flag. Flags are ordered by atom index, and by
/// issue within an atom, so the output is stable.
///
/// Runs in `O(atoms)` time.
#[must_use]
pub fn quality_flags(structure: &Structure) -> Vec<QualityFlag> {
    let mut flags = Vec::new();
    for atom in structure.data().atoms() {
        let index = atom.index();
        if let Some(occupancy) = atom.occupancy() {
            if occupancy == 0.0 {
                flags.push(QualityFlag {
                    atom: index,
                    issue: QualityIssue::ZeroOccupancy,
                });
            } else if !(0.0..=1.0).contains(&occupancy) {
                flags.push(QualityFlag {
                    atom: index,
                    issue: QualityIssue::OccupancyOutOfRange,
                });
            }
        }
        if let Some(b_factor) = atom.b_factor()
            && b_factor < 0.0
        {
            flags.push(QualityFlag {
                atom: index,
                issue: QualityIssue::NegativeBFactor,
            });
        }
    }
    flags
}

#[cfg(test)]
#[path = "quality_tests.rs"]
mod tests;
