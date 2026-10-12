//! Global distance test: the best superposition for each distance cutoff.
//!
//! GDT answers "what is the largest part of the model that can be laid within
//! `c` ångström of the reference", asked separately for every cutoff `c` and
//! averaged. The word *largest* is the whole definition: the score is a maximum
//! over rigid superpositions, taken independently per cutoff. A single
//! least-squares fit is not that maximum. It spends its freedom on the atoms
//! that deviate most, so a model with one displaced domain is fitted to a
//! compromise that leaves few atoms inside the tight cutoffs.
//!
//! # Algorithm
//!
//! This is the iterative local-superposition search of LGA and of the CASP
//! assessment software.
//!
//! 1. Seeds are contiguous fragments of the correspondence list. Fragment
//!    lengths start at the full length and halve down to four atoms; each length
//!    slides along the list in steps of a quarter of its length (and wide enough
//!    that no length has more than [`MAX_SEEDS_PER_LENGTH`] starts).
//! 2. For a seed and a cutoff, the fragment is superposed by least squares and
//!    the atoms whose distance is within that cutoff form the next fitting set.
//!    The set is refitted and the step repeated until it stops changing, for at
//!    most [`MAX_ITERATIONS`] rounds. A set smaller than three atoms is widened
//!    by half an ångström at a time until it can be fitted.
//! 3. Every superposition visited is scored against all cutoffs, and each
//!    cutoff keeps the largest count it has seen. A superposition found while
//!    chasing one cutoff therefore also counts for the others, which keeps the
//!    counts monotone in the cutoff.
//! 4. The score is the mean over cutoffs of `best count / length`.
//!
//! The search is a heuristic for a maximisation that is not solvable exactly in
//! polynomial time, so it returns a lower bound on the true maximum, as LGA
//! does. Seeds include the full-length fit, so the result is never below the
//! plain least-squares GDT.
//!
//! # Correspondence
//!
//! The `i`-th model point is paired with the `i`-th reference point and that
//! pairing is never changed: this is LGA's sequence-dependent mode. A
//! sequence-independent search (re-pairing atoms by distance, as in LGA `-4`)
//! can score higher on structures with insertions or domain swaps and is a
//! different quantity; build the pairing first with the mapping functions.
//!
//! # Cost
//!
//! At most `4 · MAX_SEEDS_PER_LENGTH · log2(n)` seeds, each refined for at most
//! [`MAX_ITERATIONS`] rounds per cutoff, and every round is one `O(n)`
//! superposition: `O(cutoffs · log n · n)` with a fixed constant, and all
//! buffers are allocated once.

use molframe_geom::{Rigid, superpose};

use crate::failure::CompareError;
use crate::numeric::{u64_to_f64, usize_to_f64};

/// The four distance cutoffs GDT-TS (total score) averages over, in ångström.
const GDT_TS_CUTOFFS: [f64; 4] = [1.0, 2.0, 4.0, 8.0];

/// The tighter cutoffs of GDT-HA (high accuracy), each half of the GDT-TS set.
const GDT_HA_CUTOFFS: [f64; 4] = [0.5, 1.0, 2.0, 4.0];

/// Refinement rounds per seed and cutoff before the search gives up on
/// convergence. Converging searches stop long before this.
const MAX_ITERATIONS: usize = 20;

/// The most start positions tried for any one seed length.
const MAX_SEEDS_PER_LENGTH: usize = 128;

/// The shortest seed fragment, in atoms.
const MIN_SEED: usize = 4;

/// How much a fitting set smaller than three atoms is widened per attempt.
const WIDEN_STEP: f64 = 0.5;

/// Computes GDT-TS: the mean over the 1/2/4/8 Å cutoffs of the largest fraction
/// of atoms that any superposition places within the cutoff.
///
/// A perfect match scores `1.0`. See the module documentation for the search.
///
/// # Errors
///
/// Returns [`CompareError::LengthMismatch`] on unequal lengths, or
/// [`CompareError::Superpose`] when no fit exists (fewer than three points, or a
/// degenerate arrangement).
pub fn gdt_ts(model: &[[f32; 3]], reference: &[[f32; 3]]) -> Result<f64, CompareError> {
    gdt_with_cutoffs(model, reference, &GDT_TS_CUTOFFS)
}

/// Computes GDT-HA, the high-accuracy variant over the tighter 0.5/1/2/4 Å set.
///
/// A perfect match scores `1.0`.
///
/// # Errors
///
/// Returns the same errors as [`gdt_ts`].
pub fn gdt_ha(model: &[[f32; 3]], reference: &[[f32; 3]]) -> Result<f64, CompareError> {
    gdt_with_cutoffs(model, reference, &GDT_HA_CUTOFFS)
}

/// Averages, over the given cutoffs, the best fraction of atoms within each.
///
/// # Errors
///
/// Returns an error for empty, non-finite, non-positive or unordered cutoffs,
/// unequal coordinate lengths, or a failed superposition.
pub fn gdt_with_cutoffs(
    model: &[[f32; 3]],
    reference: &[[f32; 3]],
    cutoffs: &[f64],
) -> Result<f64, CompareError> {
    if cutoffs.is_empty()
        || cutoffs
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
        || cutoffs.windows(2).any(|pair| pair[0] >= pair[1])
    {
        return Err(CompareError::InvalidDistanceCutoff);
    }
    if model.len() != reference.len() {
        return Err(CompareError::LengthMismatch {
            model: model.len(),
            reference: reference.len(),
        });
    }
    // The full-length fit is the first seed; failing it reports the error the
    // caller needs (too few points, degenerate arrangement) before any search.
    let whole = superpose(model, reference).map_err(CompareError::Superpose)?;
    let mut search = Search::new(model, reference, cutoffs);
    search.visit(whole.transform);
    search.run();
    let length = usize_to_f64(reference.len());
    let fractions: f64 = search
        .best
        .iter()
        .map(|count| u64_to_f64(*count) / length)
        .sum();
    Ok(fractions / usize_to_f64(cutoffs.len()))
}

/// Working state of one GDT search. Buffers are reused across every round.
struct Search<'a> {
    model: &'a [[f32; 3]],
    reference: &'a [[f32; 3]],
    cutoffs: &'a [f64],
    /// The largest count seen so far under each cutoff.
    best: Vec<u64>,
    /// Distance of every atom under the transform being scored.
    distances: Vec<f64>,
    /// Indices of the current fitting set.
    current: Vec<usize>,
    /// Indices of the next fitting set.
    next: Vec<usize>,
    /// Gathered coordinates of the fitting set.
    moving: Vec<[f32; 3]>,
    fixed: Vec<[f32; 3]>,
}

impl<'a> Search<'a> {
    fn new(model: &'a [[f32; 3]], reference: &'a [[f32; 3]], cutoffs: &'a [f64]) -> Self {
        let length = reference.len();
        Self {
            model,
            reference,
            cutoffs,
            best: vec![0; cutoffs.len()],
            distances: vec![0.0; length],
            current: Vec::with_capacity(length),
            next: Vec::with_capacity(length),
            moving: Vec::with_capacity(length),
            fixed: Vec::with_capacity(length),
        }
    }

    /// Scores one transform against every cutoff and keeps the maxima.
    fn visit(&mut self, transform: Rigid) {
        let mut within = vec![0u64; self.cutoffs.len()];
        for ((slot, &point), &target) in self
            .distances
            .iter_mut()
            .zip(self.model)
            .zip(self.reference)
        {
            let distance = molframe_geom::distance(transform.apply(point), target);
            *slot = distance;
            for (cutoff, count) in self.cutoffs.iter().zip(&mut within) {
                if distance <= *cutoff {
                    *count += 1;
                }
            }
        }
        for (best, count) in self.best.iter_mut().zip(within) {
            *best = (*best).max(count);
        }
    }

    /// Fits the atoms in `current`, or `None` when they cannot be fitted.
    fn fit_current(&mut self) -> Option<Rigid> {
        self.moving.clear();
        self.fixed.clear();
        for &index in &self.current {
            self.moving.push(self.model[index]);
            self.fixed.push(self.reference[index]);
        }
        superpose(&self.moving, &self.fixed)
            .ok()
            .map(|fit| fit.transform)
    }

    fn run(&mut self) {
        let length = self.reference.len();
        let mut fragment = length;
        loop {
            let last_start = length - fragment;
            let step = (fragment / 4)
                .max(1)
                .max(last_start.div_ceil(MAX_SEEDS_PER_LENGTH));
            let mut start = 0;
            while start <= last_start {
                self.refine_seed(start, fragment);
                start += step;
            }
            // The last window is always tried, whatever the step did.
            if !last_start.is_multiple_of(step) {
                self.refine_seed(last_start, fragment);
            }
            if fragment / 2 < MIN_SEED {
                break;
            }
            fragment /= 2;
        }
    }

    /// Refines one seed fragment once per cutoff.
    fn refine_seed(&mut self, start: usize, fragment: usize) {
        for cutoff_index in 0..self.cutoffs.len() {
            self.current.clear();
            self.current.extend(start..start + fragment);
            for _ in 0..MAX_ITERATIONS {
                let Some(transform) = self.fit_current() else {
                    break;
                };
                self.visit(transform);
                if !self.select_next(self.cutoffs[cutoff_index]) || self.next == self.current {
                    break;
                }
                std::mem::swap(&mut self.current, &mut self.next);
            }
        }
    }

    /// Fills `next` with the atoms within `cutoff` of the last scored transform,
    /// widening the cutoff until at least three atoms qualify. Returns whether a
    /// set was found.
    fn select_next(&mut self, cutoff: f64) -> bool {
        let mut limit = cutoff;
        let ceiling = self.distances.iter().copied().fold(0.0_f64, f64::max);
        loop {
            self.next.clear();
            for (index, distance) in self.distances.iter().enumerate() {
                if *distance <= limit {
                    self.next.push(index);
                }
            }
            if self.next.len() >= 3 {
                return true;
            }
            if limit >= ceiling {
                return false;
            }
            limit += WIDEN_STEP;
        }
    }
}

#[cfg(test)]
#[path = "gdt_tests.rs"]
mod tests;
