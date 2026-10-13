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
//!    that no length has more than 128 starts).
//! 2. For a seed and a cutoff, the fragment is superposed by least squares and
//!    the atoms whose distance is within that cutoff form the next fitting set.
//!    The set is refitted and the step repeated until it stops changing, for at
//!    most 20 rounds. A set smaller than three atoms is widened
//!    by half an ångström at a time until it can be fitted.
//! 3. Every superposition visited is scored against all cutoffs, and each
//!    cutoff keeps the largest count it has seen. A superposition found while
//!    chasing one cutoff therefore also counts for the others, which keeps the
//!    counts monotone in the cutoff.
//! 4. The score is the mean over cutoffs of `best count / length`.
//!
//! The search is a heuristic for a maximisation that is not solvable exactly in
//! polynomial time, so it returns a lower bound on the true maximum, as LGA
//! does. Every returned count is a real count under a real transform on all
//! atoms, so the result never exceeds the true optimum, and the full-length fit
//! is always a candidate, so it is never below the plain least-squares GDT.
//!
//! # Correspondence
//!
//! The `i`-th model point is paired with the `i`-th reference point and that
//! pairing is never changed: this is LGA's sequence-dependent mode. A
//! sequence-independent search (re-pairing atoms by distance, as in LGA `-4`)
//! can score higher on structures with insertions or domain swaps and is a
//! different quantity; build the pairing first with the mapping functions.
//!
//! # Cost and the monotone relation
//!
//! The search lives in [`crate::seed_search`]: seeds are refined on at most
//! `SAMPLE_CAP` (1024) atoms, so seed work is a fixed amount, and only the best
//! few seed transforms per cutoff are refined and counted on every atom, which is
//! `O(n)`. Seeds run in parallel over rayon workers and are reduced in seed
//! order, so the result does not depend on the worker count.
//!
//! For up to `SAMPLE_CAP` atoms the seed set is a superset of the one used by the
//! earlier serial implementation (same fragment schedule, same refinement), and
//! GDT is a maximum, so a search that visits more superpositions can only raise
//! the score: the new value is never below the old one there. Beyond that length
//! the sampled search is a different heuristic and is only guaranteed to be at
//! most the true optimum and at least the full-length fit.

use molframe_geom::superpose;

use crate::failure::CompareError;
use crate::numeric::usize_to_f64;
use crate::seed_search::{Objective, maximise, per_atom};

/// The four distance cutoffs GDT-TS (total score) averages over, in ångström.
const GDT_TS_CUTOFFS: [f64; 4] = [1.0, 2.0, 4.0, 8.0];

/// The tighter cutoffs of GDT-HA (high accuracy), each half of the GDT-TS set.
const GDT_HA_CUTOFFS: [f64; 4] = [0.5, 1.0, 2.0, 4.0];

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
    // The full-length fit is the first candidate; failing it reports the error
    // the caller needs (too few points, degenerate arrangement) before any
    // search.
    let whole = superpose(model, reference).map_err(CompareError::Superpose)?;
    let objective = Objective::gdt(cutoffs);
    let outcome = maximise(model, reference, &objective, whole.transform);
    let fractions: f64 = outcome
        .scores
        .iter()
        .map(|count| per_atom(*count, reference.len()))
        .sum();
    Ok(fractions / usize_to_f64(cutoffs.len()))
}

#[cfg(test)]
#[path = "gdt_tests.rs"]
mod tests;
