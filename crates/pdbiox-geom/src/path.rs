//! Torsions along an ordered atom path.

use crate::dihedral;

/// Computes consecutive torsions for an atom path, capped at `maximum`.
///
/// A missing or degenerate quartet yields `None` only for the dependent
/// torsion. The function never bridges across an absent atom.
#[must_use]
pub fn path_torsions(atoms: &[Option<[f32; 3]>], maximum: usize) -> Vec<Option<f64>> {
    atoms
        .windows(4)
        .take(maximum)
        .map(|quartet| match quartet {
            [Some(a), Some(b), Some(c), Some(d)] => dihedral(*a, *b, *c, *d),
            _ => None,
        })
        .collect()
}

#[cfg(test)]
#[path = "path_tests.rs"]
mod tests;
