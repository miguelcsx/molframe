//! Scores measured after fitting one structure onto the other.
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
use crate::numeric::{u64_to_f64, usize_to_f64};

/// The four distance cutoffs GDT-TS (total score) averages over, in ångström.
const GDT_TS_CUTOFFS: [f64; 4] = [1.0, 2.0, 4.0, 8.0];

/// The tighter cutoffs of GDT-HA (high accuracy), each half of the GDT-TS set.
const GDT_HA_CUTOFFS: [f64; 4] = [0.5, 1.0, 2.0, 4.0];

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
    let fit = fit(model, reference)?;
    let length = reference.len();
    let scale = tm_scale(length);
    let mut total = 0.0f64;
    for (moved, &target) in model
        .iter()
        .map(|&point| fit.transform.apply(point))
        .zip(reference)
    {
        let distance = molframe_geom::distance(moved, target);
        let ratio = distance / scale;
        total += 1.0 / (1.0 + ratio * ratio);
    }
    Ok(total / usize_to_f64(length))
}

/// Computes GDT-TS: the mean over the 1/2/4/8 Å cutoffs of the fraction within.
///
/// A perfect match scores `1.0`. The model is fitted onto the reference first.
///
/// # Errors
///
/// Returns the same errors as [`tm_score`].
pub fn gdt_ts(model: &[[f32; 3]], reference: &[[f32; 3]]) -> Result<f64, CompareError> {
    gdt_with_cutoffs(model, reference, &GDT_TS_CUTOFFS)
}

/// Computes GDT-HA, the high-accuracy variant over the tighter 0.5/1/2/4 Å set.
///
/// A perfect match scores `1.0`. The model is fitted onto the reference first.
///
/// # Errors
///
/// Returns the same errors as [`tm_score`].
pub fn gdt_ha(model: &[[f32; 3]], reference: &[[f32; 3]]) -> Result<f64, CompareError> {
    gdt_with_cutoffs(model, reference, &GDT_HA_CUTOFFS)
}

/// Averages, over the given cutoffs, the fraction of fitted atoms within each.
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
    let fit = fit(model, reference)?;
    let length = reference.len();
    let mut within = vec![0u64; cutoffs.len()];
    for (moved, &target) in model
        .iter()
        .map(|&point| fit.transform.apply(point))
        .zip(reference)
    {
        let distance = molframe_geom::distance(moved, target);
        for (cutoff, count) in cutoffs.iter().zip(&mut within) {
            if distance <= *cutoff {
                *count += 1;
            }
        }
    }
    let fractions: f64 = within
        .iter()
        .map(|count| u64_to_f64(*count) / usize_to_f64(length))
        .sum();
    Ok(fractions / usize_to_f64(cutoffs.len()))
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
    if model.len() != reference.len() || model.len() != weights.len() {
        return Err(CompareError::LengthMismatch {
            model: model.len(),
            reference: reference.len(),
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
