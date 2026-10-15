//! The one definition of which alternate-conformation labels coexist.
//!
//! Anything that pairs atoms (clash detection, contact search, geometry
//! checks, bonding) must agree on this rule, or atoms from different
//! conformers get paired as if they were present together.

use crate::symbol::AltId;

/// Whether two atoms with these alternate labels can be present together.
///
/// A blank label marks an atom present in every conformation, so it is
/// compatible with any label, blank or not. Two non-blank labels are
/// compatible only when they are the same label; different non-blank labels
/// are different conformers of the same site and never coexist.
///
/// An atom without any recorded label should be passed as [`AltId::BLANK`].
///
/// Runs in `O(1)` time and allocates no memory.
#[must_use]
#[inline]
pub fn altloc_compatible(left: AltId, right: AltId) -> bool {
    left.is_blank() || right.is_blank() || left == right
}

/// The index pairs `(i, j)` with `i < j` whose labels are compatible.
///
/// Pairs come out in ascending `(i, j)` order, so results are deterministic.
/// Runs in `O(n²)` time over the labels and allocates nothing; use it to
/// filter candidate pairs when `n` is small (one residue, one contact cell),
/// and [`altloc_compatible`] directly when pairs come from a spatial search.
pub fn compatible_pairs(labels: &[AltId]) -> impl Iterator<Item = (usize, usize)> + '_ {
    (0..labels.len()).flat_map(move |first| {
        ((first + 1)..labels.len())
            .filter(move |second| altloc_compatible(labels[first], labels[*second]))
            .map(move |second| (first, second))
    })
}

#[cfg(test)]
#[path = "compat_tests.rs"]
mod tests;
