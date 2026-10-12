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
//! producing a reflection. Cost is one superposition plus one pass over the
//! atoms.

use molframe_geom::{Superposition, superpose};

use crate::failure::CompareError;
use crate::numeric::usize_to_f64;

pub use crate::gdt::{gdt_ha, gdt_ts, gdt_with_cutoffs};

/// Computes the TM-score of a model against a reference.
///
/// The score runs from just above zero for unrelated structures to `1.0` for
/// identical ones, normalised by the reference length so it does not depend on
/// how long the structures are. The model is fitted onto the reference first.
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
    let mut best_score = tm_sum(fit(model, reference)?.transform, model, reference, scale);

    // Reweighted fits. The gradient of `sum 1/(1 + (d/d0)^2)` with respect to a
    // point's displacement points along `-d` with magnitude `1/(1+(d/d0)^2)^2`,
    // so a weighted least-squares fit with those weights moves the fit toward
    // the TM optimum. Reweighting changes the objective, so the loop keeps the
    // best score rather than trusting the last fit, and stops as soon as a step
    // fails to improve it.
    //
    // This is a local search, while the published score is a maximum over
    // superpositions. TM-align reaches that maximum with a residue-level dynamic
    // program that can re-pair residues as well as move them, which a
    // correspondence-preserving kernel cannot do, so on a structure whose parts
    // have swung past each other it reports a little more than this does — about
    // 0.04 for a 1UBQ half-and-half hinge, and about 0.015 for a displaced tail.
    // The value here is the best of the superpositions this kernel can reach, it
    // is never below the plain least-squares fit, and it is exact for a rigid
    // displacement.
    let mut candidate = fit(model, reference)?.transform;
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
    Ok(best_score / usize_to_f64(length))
}

/// How many reweighted fits the TM search takes before it stops improving.
///
/// Each step is one quaternion decomposition, so the whole search costs a small
/// constant number of them; the loop leaves as soon as a step stops improving.
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
