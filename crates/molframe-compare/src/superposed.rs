//! Scores measured after fitting one structure onto the other.
//!
//! GDT lives in its own module because it is a search over superpositions, not
//! one fit; it is re-exported here.
//!
//! TM-score and GDT both quantify how close corresponding atoms come once the
//! model is optimally superposed on the reference. They differ in shape: GDT is
//! the mean fraction of atoms under a few fixed distance cutoffs, while TM-score
//! weights every atom smoothly by a length-dependent scale, so that adding
//! disordered tails cannot drag a good core score down.
//!
//! Both use the shared quaternion superposition, which fits without ever
//! producing a reflection, and both maximise over superpositions with the
//! fragment-seeded search in `seed_search`. Its cost is bounded (seeds refine on
//! at most 1024 atoms) plus a few `O(n)` passes, and it runs in parallel
//! deterministically.

use molframe_geom::{Superposition, superpose};
use rayon::prelude::*;

use crate::failure::CompareError;
use crate::numeric::usize_to_f64;
use crate::seed_search::{Objective, maximise};

pub use crate::gdt::{gdt_ha, gdt_ts, gdt_with_cutoffs};

/// Computes the TM-score of a model against a reference.
///
/// The score runs from just above zero for unrelated structures to `1.0` for
/// identical ones, normalised by the reference length so it does not depend on
/// how long the structures are. The model is fitted onto the reference first.
///
/// The score is the maximum over rigid superpositions, searched from several
/// fragment-seeded starting fits (see the module documentation), so a
/// multi-domain or hinged pair is scored by its best-fitting domain arrangement
/// rather than by one global least-squares compromise. The correspondence is
/// the given one: residue `i` of the model answers to residue `i` of the
/// reference and is never re-paired. TM-align additionally re-pairs residues by
/// dynamic programming, so on pairs where that matters it can report more; this
/// function, given the same correspondence, approaches the TM-score program's
/// optimum. The value is always the score of an actual superposition, so it
/// never exceeds the true maximum for the correspondence.
///
/// # Errors
///
/// Returns [`CompareError::LengthMismatch`] on unequal lengths, or
/// [`CompareError::Superpose`] when no fit exists (fewer than three points, or a
/// degenerate arrangement).
pub fn tm_score(model: &[[f32; 3]], reference: &[[f32; 3]]) -> Result<f64, CompareError> {
    if model.len() != reference.len() {
        return Err(CompareError::LengthMismatch {
            model: model.len(),
            reference: reference.len(),
        });
    }
    let length = reference.len();
    let scale = tm_scale(length);
    let whole = fit(model, reference)?.transform;
    let mut best_score = tm_sum(whole, model, reference, scale);

    // Multi-seed search, as in TM-score and TM-align: superpositions seeded from
    // fragments of several lengths, each refined by repeatedly keeping the atoms
    // near their targets and refitting, scored by the TM sum. See
    // `crate::seed_search` for the seed schedule, the bound on its cost and why
    // the result does not depend on the number of workers.
    let outcome = maximise(model, reference, &Objective::tm(scale), whole);
    best_score = best_score.max(outcome.scores[0]);

    // Reweighted polish of the leading transforms. The gradient of
    // `sum 1/(1 + (d/d0)^2)` with respect to a point's displacement points along
    // `-d` with magnitude `1/(1+(d/d0)^2)^2`, so a weighted least-squares fit
    // with those weights moves a fit toward a local TM optimum. Each step is kept
    // only if it improves the score, so this never lowers the result.
    let starts: Vec<molframe_geom::Rigid> = std::iter::once(whole)
        .chain(outcome.leaders.iter().copied())
        .collect();
    let polished: Vec<f64> = starts
        .par_iter()
        .map(|&start| reweighted_score(start, model, reference, scale))
        .collect();
    for score in polished {
        best_score = best_score.max(score);
    }
    Ok(best_score / usize_to_f64(length))
}

/// The best TM sum reached by reweighted fits that start from `start`.
fn reweighted_score(
    start: molframe_geom::Rigid,
    model: &[[f32; 3]],
    reference: &[[f32; 3]],
    scale: f64,
) -> f64 {
    let mut candidate = start;
    let mut best_score = tm_sum(candidate, model, reference, scale);
    for _ in 0..TM_REWEIGHT_STEPS {
        let Some(weights) = tm_weights(candidate, model, reference, scale) else {
            break;
        };
        let Ok(step) = molframe_geom::superpose_weighted(model, reference, &weights) else {
            break;
        };
        let score = tm_sum(step.transform, model, reference, scale);
        if score <= best_score {
            break;
        }
        candidate = step.transform;
        best_score = score;
    }
    best_score
}

/// How many reweighted fits the polish takes before it stops improving.
const TM_REWEIGHT_STEPS: usize = 32;

/// The TM sum for one transform, before normalisation.
fn tm_sum(
    transform: molframe_geom::Rigid,
    model: &[[f32; 3]],
    reference: &[[f32; 3]],
    scale: f64,
) -> f64 {
    let mut total = 0.0f64;
    for (moved, &target) in model
        .iter()
        .map(|&point| transform.apply(point))
        .zip(reference)
    {
        let ratio = molframe_geom::distance(moved, target) / scale;
        total += 1.0 / (1.0 + ratio * ratio);
    }
    total
}

/// The per-point weights of one reweighted fit, or `None` when they carry no
/// information.
///
/// A weight is the squared TM kernel, so a point already close to its target
/// pulls hardest. All-zero weights leave the fit undetermined, which is the
/// signal to stop rather than to fit noise.
fn tm_weights(
    transform: molframe_geom::Rigid,
    model: &[[f32; 3]],
    reference: &[[f32; 3]],
    scale: f64,
) -> Option<Vec<f64>> {
    let mut weights = Vec::with_capacity(model.len());
    let mut total = 0.0f64;
    for (moved, &target) in model
        .iter()
        .map(|&point| transform.apply(point))
        .zip(reference)
    {
        let ratio = molframe_geom::distance(moved, target) / scale;
        let kernel = 1.0 / (1.0 + ratio * ratio);
        let weight = kernel * kernel;
        total += weight;
        weights.push(weight);
    }
    (total > 0.0 && total.is_finite()).then_some(weights)
}

/// Computes weighted RMSD for coordinates already in the chosen frame.
///
/// Weights must be finite and non-negative, with at least one positive value.
/// Alignment remains a separate explicit workflow stage.
///
/// # Errors
///
/// Returns a length mismatch or invalid-score-input error.
pub fn weighted_rmsd(
    model: &[[f32; 3]],
    reference: &[[f32; 3]],
    weights: &[f64],
) -> Result<f64, CompareError> {
    if model.len() != reference.len() {
        return Err(CompareError::LengthMismatch {
            model: model.len(),
            reference: reference.len(),
        });
    }
    if model.len() != weights.len() {
        return Err(CompareError::WeightLengthMismatch {
            points: model.len(),
            weights: weights.len(),
        });
    }
    let mut weighted_squared = 0.0;
    let mut total_weight = 0.0;
    for ((model, reference), weight) in model.iter().zip(reference).zip(weights) {
        if !weight.is_finite() || *weight < 0.0 {
            return Err(CompareError::InvalidScoreInput);
        }
        weighted_squared += weight * molframe_geom::distance_squared(*model, *reference);
        total_weight += weight;
    }
    if total_weight <= 0.0 || !total_weight.is_finite() {
        return Err(CompareError::InvalidScoreInput);
    }
    Ok((weighted_squared / total_weight).sqrt())
}

/// Superposes the model onto the reference, translating the error type.
fn fit(model: &[[f32; 3]], reference: &[[f32; 3]]) -> Result<Superposition, CompareError> {
    if model.len() != reference.len() {
        return Err(CompareError::LengthMismatch {
            model: model.len(),
            reference: reference.len(),
        });
    }
    superpose(model, reference).map_err(CompareError::Superpose)
}

/// The length-dependent distance scale d0 used by the TM-score.
fn tm_scale(length: usize) -> f64 {
    if length > 21 {
        let cube_root = (usize_to_f64(length) - 15.0).cbrt();
        (1.24 * cube_root - 1.8).max(0.5)
    } else {
        0.5
    }
}

#[cfg(test)]
#[path = "superposed_tests.rs"]
mod tests;
