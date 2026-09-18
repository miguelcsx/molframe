//! How much two contact maps overlap.
//!
//! Given the contacting residue pairs of two structures, the Jaccard index —
//! the shared pairs over the pairs in either — says how similar their contact
//! patterns are without any superposition. Each pair is read unordered, so it
//! does not matter which residue a caller listed first.

use std::collections::BTreeSet;

use crate::numeric::usize_to_f64;

/// The overlap between two sets of contacts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactSimilarity {
    /// Pairs present in both maps.
    pub shared: usize,
    /// Pairs present in either map.
    pub union: usize,
    /// `shared / union`, or `1.0` when both maps are empty.
    pub jaccard: f64,
}

/// Compares two contact maps given as lists of residue pairs.
///
/// Pairs are treated as unordered and deduplicated, so a repeated or reversed
/// pair does not change the result. Two empty maps are perfectly similar.
///
/// Runs in `O(n log n)` in the total number of pairs.
#[must_use]
pub fn contact_map_similarity(first: &[(u32, u32)], second: &[(u32, u32)]) -> ContactSimilarity {
    let left = normalise(first);
    let right = normalise(second);
    let shared = left.intersection(&right).count();
    let union = left.union(&right).count();
    let jaccard = if union == 0 {
        1.0
    } else {
        usize_to_f64(shared) / usize_to_f64(union)
    };
    ContactSimilarity {
        shared,
        union,
        jaccard,
    }
}

/// Collects pairs into an unordered, deduplicated set.
fn normalise(pairs: &[(u32, u32)]) -> BTreeSet<(u32, u32)> {
    pairs
        .iter()
        .map(|&(a, b)| if a <= b { (a, b) } else { (b, a) })
        .collect()
}

#[cfg(test)]
#[path = "similarity_tests.rs"]
mod tests;
