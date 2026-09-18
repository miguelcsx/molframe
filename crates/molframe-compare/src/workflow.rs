//! Explicit mapping, alignment, measurement and verdict stages.

use crate::CompareError;
use crate::numeric::usize_to_f64;
use molframe_geom::{Rigid, distance, superpose};
use std::collections::BTreeSet;

/// One validated point correspondence.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct PointMatch {
    /// Point in the reference coordinate set.
    pub reference: usize,
    /// Corresponding point in the model coordinate set.
    pub model: usize,
}

/// Stable, one-to-one correspondence established before fitting or scoring.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PointMapping(Box<[PointMatch]>);

impl PointMapping {
    /// Validates bounds and one-to-one use, then preserves declaration order.
    ///
    /// # Errors
    ///
    /// Returns a mapping error for an empty, out-of-range or repeated point.
    pub fn new(
        matches: impl IntoIterator<Item = PointMatch>,
        reference_len: usize,
        model_len: usize,
    ) -> Result<Self, CompareError> {
        let matches: Vec<_> = matches.into_iter().collect();
        let mut references = BTreeSet::new();
        let mut models = BTreeSet::new();
        if matches.is_empty()
            || matches.iter().any(|mapping| {
                mapping.reference >= reference_len
                    || mapping.model >= model_len
                    || !references.insert(mapping.reference)
                    || !models.insert(mapping.model)
            })
        {
            return Err(CompareError::InvalidMapping);
        }
        Ok(Self(matches.into_boxed_slice()))
    }

    /// Correspondences in the caller's stable order.
    #[must_use]
    pub fn matches(&self) -> &[PointMatch] {
        &self.0
    }
}

/// Alignment explicitly chosen for a mapped comparison.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ComparisonAlignment {
    /// Measure in the input coordinate frame.
    NotRequired,
    /// Apply this rigid transform to model coordinates first.
    Rigid(Rigid),
}

/// Fits the mapped model points to the mapped reference points.
///
/// # Errors
///
/// Returns a superposition error when the mapping cannot fix a rigid pose.
pub fn align_mapping(
    reference: &[[f32; 3]],
    model: &[[f32; 3]],
    mapping: &PointMapping,
) -> Result<ComparisonAlignment, CompareError> {
    let reference: Vec<_> = mapping
        .matches()
        .iter()
        .map(|pair| reference[pair.reference])
        .collect();
    let model: Vec<_> = mapping
        .matches()
        .iter()
        .map(|pair| model[pair.model])
        .collect();
    superpose(&model, &reference)
        .map(|fit| ComparisonAlignment::Rigid(fit.transform))
        .map_err(CompareError::Superpose)
}

/// Raw per-pair distances and their RMS aggregate.
#[derive(Clone, Debug, PartialEq)]
pub struct DistanceMeasurement {
    /// Distance for each mapping row, in ångström.
    pub distances: Box<[f64]>,
    /// Root-mean-square mapped distance, in ångström.
    pub rmsd: f64,
}

/// Measures a fixed mapping after applying the separately selected alignment.
#[must_use]
pub fn measure_mapping(
    reference: &[[f32; 3]],
    model: &[[f32; 3]],
    mapping: &PointMapping,
    alignment: ComparisonAlignment,
) -> DistanceMeasurement {
    let distances: Box<[f64]> = mapping
        .matches()
        .iter()
        .map(|pair| {
            let model = match alignment {
                ComparisonAlignment::NotRequired => model[pair.model],
                ComparisonAlignment::Rigid(transform) => transform.apply(model[pair.model]),
            };
            distance(model, reference[pair.reference])
        })
        .collect::<Vec<_>>()
        .into_boxed_slice();
    let rmsd = (distances.iter().map(|value| value * value).sum::<f64>()
        / usize_to_f64(distances.len()))
    .sqrt();
    DistanceMeasurement { distances, rmsd }
}

/// Threshold decision kept separate from the raw measurement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ComparisonVerdict {
    /// Versioned-profile or caller-provided threshold in ångström.
    pub maximum_rmsd: f64,
    /// Whether the raw RMSD is at or below the threshold.
    pub passed: bool,
}

/// Applies a threshold without discarding the underlying measurement.
#[must_use]
pub fn decide_rmsd(measurement: &DistanceMeasurement, maximum_rmsd: f64) -> ComparisonVerdict {
    ComparisonVerdict {
        maximum_rmsd,
        passed: maximum_rmsd.is_finite() && measurement.rmsd <= maximum_rmsd,
    }
}

#[cfg(test)]
#[path = "workflow_tests.rs"]
mod tests;
