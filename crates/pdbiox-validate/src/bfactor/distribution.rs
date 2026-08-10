//! Distribution summaries and explicit standard-score outliers.

use pdbiox_core::index::AtomIndex;
use pdbiox_core::{AtomSelection, Structure};

use crate::numeric::usize_to_f64;

/// Distribution summary and atom-level outliers.
#[derive(Clone, Debug, PartialEq)]
pub struct BFactorDistribution {
    /// Selected atoms intended for assessment.
    pub intended: usize,
    /// Selected atoms with recorded B factors.
    pub assessed: usize,
    /// Arithmetic mean.
    pub mean: f64,
    /// Population variance.
    pub variance: f64,
    /// Population standard deviation.
    pub standard_deviation: f64,
    /// Minimum recorded value.
    pub minimum: f64,
    /// Median recorded value.
    pub median: f64,
    /// Maximum recorded value.
    pub maximum: f64,
    /// Values beyond the caller's explicit absolute z-score threshold.
    pub outliers: Vec<BFactorOutlier>,
}

/// One B-factor distribution outlier.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BFactorOutlier {
    /// Atom index.
    pub atom: AtomIndex,
    /// Recorded B factor.
    pub value: f64,
    /// Population standard score.
    pub z_score: f64,
}

/// Invalid B-factor analysis controls or absent observations.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum BFactorError {
    /// Z-score threshold must be finite and non-negative.
    #[error("B-factor z-score threshold must be finite and non-negative")]
    InvalidThreshold,
    /// No selected atom has a recorded B factor.
    #[error("B-factor distribution requires at least one recorded selected value")]
    NoObservations,
    /// TLS matrices, tolerance or coordinates are invalid.
    #[error("TLS inputs must be finite with symmetric T/L matrices and non-negative tolerance")]
    InvalidTls,
}

/// Summarizes selected B factors and reports explicit z-score outliers.
///
/// # Errors
///
/// Returns invalid-threshold or no-observation errors.
pub fn b_factor_distribution(
    structure: &Structure,
    selection: &AtomSelection,
    outlier_standard_deviations: f64,
) -> Result<BFactorDistribution, BFactorError> {
    if !outlier_standard_deviations.is_finite() || outlier_standard_deviations < 0.0 {
        return Err(BFactorError::InvalidThreshold);
    }
    let values: Vec<_> = selection
        .iter()
        .filter_map(|index| {
            let atom = structure.data().atom(AtomIndex::new(index))?;
            atom.b_factor()
                .map(|value| (atom.index(), f64::from(value)))
        })
        .collect();
    if values.is_empty() {
        return Err(BFactorError::NoObservations);
    }
    let count = usize_to_f64(values.len());
    let mean = values.iter().map(|(_, value)| value).sum::<f64>() / count;
    let variance = values
        .iter()
        .map(|(_, value)| (value - mean).powi(2))
        .sum::<f64>()
        / count;
    let standard_deviation = variance.sqrt();
    let mut sorted: Vec<f64> = values.iter().map(|(_, value)| *value).collect();
    sorted.sort_by(f64::total_cmp);
    let median = if sorted.len().is_multiple_of(2) {
        f64::midpoint(sorted[sorted.len() / 2 - 1], sorted[sorted.len() / 2])
    } else {
        sorted[sorted.len() / 2]
    };
    let outliers = if standard_deviation > 0.0 {
        values
            .iter()
            .filter_map(|(atom, value)| {
                let z_score = (value - mean) / standard_deviation;
                (z_score.abs() > outlier_standard_deviations).then_some(BFactorOutlier {
                    atom: *atom,
                    value: *value,
                    z_score,
                })
            })
            .collect()
    } else {
        Vec::new()
    };
    Ok(BFactorDistribution {
        intended: selection.iter().count(),
        assessed: values.len(),
        mean,
        variance,
        standard_deviation,
        minimum: sorted[0],
        median,
        maximum: sorted[sorted.len() - 1],
        outliers,
    })
}
