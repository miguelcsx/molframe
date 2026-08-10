//! Counting the fixed-length words of a sequence.
//!
//! A k-mer table is the multiset of every length-`k` window of a sequence. It is
//! the starting point for seed-and-extend search and for cheap composition
//! comparisons, and it is deterministic: the words come back in sorted order, so
//! two runs over the same sequence agree exactly.
//!
//! Cost is `O(n · k)` time and `O(distinct k-mers · k)` space.

use std::collections::BTreeMap;

/// Counts every length-`k` window of `sequence`, sorted by the word.
///
/// A `k` of zero, or one longer than the sequence, yields an empty table since
/// there are no windows to count.
///
/// # Examples
///
/// ```
/// use pdbiox_seq::kmer_counts;
///
/// let counts = kmer_counts(b"AAAA", 2);
/// assert_eq!(counts, vec![(b"AA".to_vec(), 3)]);
/// ```
#[must_use]
pub fn kmer_counts(sequence: &[u8], k: usize) -> Vec<(Vec<u8>, u32)> {
    if k == 0 || sequence.len() < k {
        return Vec::new();
    }
    let mut counts: BTreeMap<&[u8], u32> = BTreeMap::new();
    for window in sequence.windows(k) {
        *counts.entry(window).or_insert(0) += 1;
    }
    counts
        .into_iter()
        .map(|(word, count)| (word.to_vec(), count))
        .collect()
}

/// Selects the minimizers of a sequence: the smallest k-mer in each window.
///
/// A minimizer scheme slides a window of `window` consecutive k-mers and keeps
/// the lexicographically smallest one — ties broken toward the leftmost — as a
/// compact, order-stable sample of the sequence. Consecutive windows that pick
/// the same position are reported once. Two sequences that share a stretch share
/// its minimizers, which is what makes them useful as anchors.
///
/// Returns the chosen `(position, k-mer)` pairs in order. A `k` or `window` of
/// zero, or a sequence too short to fill one window, yields an empty result.
///
/// Runs in `O(kmers · window · k)` time.
#[must_use]
pub fn minimizers(sequence: &[u8], k: usize, window: usize) -> Vec<(usize, Vec<u8>)> {
    if k == 0 || window == 0 || sequence.len() < k {
        return Vec::new();
    }
    let kmer_count = sequence.len() - k + 1;
    if kmer_count < window {
        return Vec::new();
    }

    let mut chosen = Vec::new();
    let mut last_position: Option<usize> = None;
    for start in 0..=(kmer_count - window) {
        let mut best = start;
        for offset in 1..window {
            let candidate = start + offset;
            if sequence[candidate..candidate + k] < sequence[best..best + k] {
                best = candidate;
            }
        }
        if last_position != Some(best) {
            chosen.push((best, sequence[best..best + k].to_vec()));
            last_position = Some(best);
        }
    }
    chosen
}

#[cfg(test)]
#[path = "counts_tests.rs"]
mod tests;
