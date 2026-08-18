//! The local distance difference test.
//!
//! lDDT asks a superposition-free question: of all the interatomic distances a
//! reference structure has within a neighbourhood, how many does the model
//! reproduce? A distance is scored against four tolerances, and an atom pair's
//! score is the fraction of tolerances it satisfies; the structure's lDDT is the
//! mean of that over every reference distance short enough to be local.
//!
//! Because it never fits one structure onto the other, no alignment can flatter
//! the score, and a rigid-body motion of part of the structure is penalised
//! exactly where it breaks the local distances. Cost is `O(n²)` in the number of
//! points, which is what comparing all short-range distances requires.

use crate::failure::CompareError;
use crate::numeric::{u64_to_f64, usize_to_f64};

/// Result policy when no reference pair belongs to the requested local domain.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EmptyLddtPolicy {
    /// Standard lDDT convention: an empty comparison is vacuously perfect.
    Perfect,
    /// Treat the absence of eligible pairs as an error.
    Error,
}

/// Explicit controls for the general local-distance preservation kernel.
#[derive(Clone, Debug, PartialEq)]
pub struct LddtOptions {
    /// Largest reference distance included in the local domain, in ångström.
    pub inclusion_radius: f64,
    /// Reference distances at or below this value are excluded.
    pub minimum_reference_distance: f64,
    /// Absolute-error tolerances averaged for every eligible pair.
    pub tolerances: Box<[f64]>,
    /// Behavior when the requested domain has no eligible pairs.
    pub empty_policy: EmptyLddtPolicy,
}

impl LddtOptions {
    /// Standard global lDDT profile with 0.5/1/2/4 Å tolerances.
    #[must_use]
    pub fn standard(inclusion_radius: f64) -> Self {
        Self {
            inclusion_radius,
            minimum_reference_distance: 0.0,
            tolerances: Box::new([0.5, 1.0, 2.0, 4.0]),
            empty_policy: EmptyLddtPolicy::Perfect,
        }
    }
}

/// Computes the global lDDT of a model against a reference.
///
/// Only reference distances at or below `inclusion_radius` count, which is what
/// makes the test local. A model that reproduces every such distance scores
/// `1.0`; one that preserves none scores `0.0`. When the reference has no local
/// distances at all, the score is `1.0`, since there is nothing to get wrong.
///
/// # Errors
///
/// Returns [`CompareError::LengthMismatch`] when the two sets differ in length.
pub fn lddt(
    model: &[[f32; 3]],
    reference: &[[f32; 3]],
    inclusion_radius: f32,
) -> Result<f64, CompareError> {
    lddt_with_options(
        model,
        reference,
        &LddtOptions::standard(f64::from(inclusion_radius)),
    )
}

/// Computes local-distance preservation using caller-declared domain and tolerances.
///
/// # Errors
///
/// Returns an error for unequal lengths, invalid coordinates/options, or an
/// empty domain when [`EmptyLddtPolicy::Error`] is selected.
pub fn lddt_with_options(
    model: &[[f32; 3]],
    reference: &[[f32; 3]],
    options: &LddtOptions,
) -> Result<f64, CompareError> {
    if model.len() != reference.len() {
        return Err(CompareError::LengthMismatch {
            model: model.len(),
            reference: reference.len(),
        });
    }
    validate_options(model, reference, options)?;

    let mut considered = 0u64;
    let mut preserved = 0.0f64;
    let minimum_squared = options.minimum_reference_distance * options.minimum_reference_distance;
    let inclusion_squared = options.inclusion_radius * options.inclusion_radius;
    for i in 0..reference.len() {
        for j in (i + 1)..reference.len() {
            let reference_squared = squared_distance(reference[i], reference[j]);
            if reference_squared <= minimum_squared || reference_squared > inclusion_squared {
                continue;
            }
            considered += 1;
            let reference_distance = reference_squared.sqrt();
            let model_distance = pdbiox_geom::distance(model[i], model[j]);
            let error = (model_distance - reference_distance).abs();
            let within = options.tolerances.len()
                - options
                    .tolerances
                    .partition_point(|tolerance| *tolerance <= error);
            preserved += usize_to_f64(within) / usize_to_f64(options.tolerances.len());
        }
    }

    if considered == 0 {
        return match options.empty_policy {
            EmptyLddtPolicy::Perfect => Ok(1.0),
            EmptyLddtPolicy::Error => Err(CompareError::NoComparablePairs),
        };
    }
    Ok(preserved / u64_to_f64(considered))
}

#[inline]
fn squared_distance(left: [f32; 3], right: [f32; 3]) -> f64 {
    let dx = f64::from(left[0]) - f64::from(right[0]);
    let dy = f64::from(left[1]) - f64::from(right[1]);
    let dz = f64::from(left[2]) - f64::from(right[2]);
    dx.mul_add(dx, dy.mul_add(dy, dz * dz))
}

fn validate_options(
    model: &[[f32; 3]],
    reference: &[[f32; 3]],
    options: &LddtOptions,
) -> Result<(), CompareError> {
    let valid_tolerances = !options.tolerances.is_empty()
        && options
            .tolerances
            .iter()
            .all(|value| value.is_finite() && *value > 0.0)
        && options.tolerances.windows(2).all(|pair| pair[0] < pair[1]);
    let valid_domain = options.inclusion_radius.is_finite()
        && options.inclusion_radius > 0.0
        && options.minimum_reference_distance.is_finite()
        && options.minimum_reference_distance >= 0.0
        && options.minimum_reference_distance < options.inclusion_radius;
    let finite_coordinates = model
        .iter()
        .chain(reference)
        .flatten()
        .all(|value| value.is_finite());
    if valid_tolerances && valid_domain && finite_coordinates {
        Ok(())
    } else {
        Err(CompareError::InvalidScoreInput)
    }
}

#[cfg(test)]
#[path = "lddt_tests.rs"]
mod tests;
