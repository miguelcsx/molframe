//! Seed-and-extend search.
//!
//! Aligning a short query against a long reference with full dynamic programming
//! is wasteful when almost none of the reference is relevant. Seed-and-extend
//! first checks for a shared exact word of length `k` — a seed — and only when
//! one exists does it pay for a local alignment. The seed is a fast filter; the
//! extension is the exact Smith–Waterman result, so a hit is never worse than
//! doing the alignment directly, and a miss is rejected in linear time.

use std::collections::BTreeSet;

use crate::align::{AlignError, Alignment, local};
use crate::scoring::Scoring;

/// Searches for a local alignment of `query` in `reference`, gated by a `k`-mer.
///
/// Returns `None` when the two share no exact word of length `k` — there is
/// nothing to extend — and otherwise the best local alignment between them. A
/// `k` of zero, or longer than either sequence, cannot seed and returns `None`.
///
/// Runs in `O(reference · k)` to seed, plus a local alignment on a hit.
///
/// # Errors
///
/// Returns [`AlignError::NumericOverflow`] when the exact extension alignment
/// exceeds its supported numeric domain.
pub fn seed_and_extend(
    query: &[u8],
    reference: &[u8],
    k: usize,
    scoring: Scoring,
) -> Result<Option<Alignment>, AlignError> {
    if k == 0 || query.len() < k || reference.len() < k {
        return Ok(None);
    }
    let seeds: BTreeSet<&[u8]> = reference.windows(k).collect();
    if !query.windows(k).any(|word| seeds.contains(word)) {
        return Ok(None);
    }
    local(query, reference, scoring).map(Some)
}

#[cfg(test)]
#[path = "search_tests.rs"]
mod tests;
