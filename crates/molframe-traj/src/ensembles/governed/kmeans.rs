//! Governed k-means wrappers over owned and borrowed numerical observations.

use super::ensemble::{GovernedEnsembleError, float, governed, integer};
use crate::{KMeans, KMeansOptions, kmeans, kmeans_view};
use molframe_core::contract::{Analysis, AnalysisPolicy, ParameterValue};

/// Governed deterministic k-means.
///
/// # Errors
///
/// Returns the typed raw-kernel failure or coverage overflow.
pub fn analyse_kmeans(
    observations: &[Vec<f64>],
    options: KMeansOptions<'_>,
    policy: &AnalysisPolicy,
) -> Result<Analysis<KMeans>, GovernedEnsembleError> {
    kmeans_analysis(
        kmeans(observations, options)?,
        observations.len(),
        options,
        policy,
    )
}

/// Governed deterministic k-means over borrowed row-major observations.
///
/// The input stays borrowed for the full native execution.  This is the
/// zero-copy counterpart to [`analyse_kmeans`] for contiguous numerical data.
///
/// # Errors
///
/// Returns the typed raw-kernel failure or coverage overflow.
pub fn analyse_kmeans_view(
    values: &[f64],
    rows: usize,
    columns: usize,
    options: KMeansOptions<'_>,
    policy: &AnalysisPolicy,
) -> Result<Analysis<KMeans>, GovernedEnsembleError> {
    kmeans_analysis(
        kmeans_view(values, rows, columns, options)?,
        rows,
        options,
        policy,
    )
}

fn kmeans_analysis(
    value: KMeans,
    observations: usize,
    options: KMeansOptions<'_>,
    policy: &AnalysisPolicy,
) -> Result<Analysis<KMeans>, GovernedEnsembleError> {
    let mut result = governed(value, observations, policy, "kmeans")?;
    result.provenance = result
        .provenance
        .with_parameter(
            "initial_centres",
            ParameterValue::Text(format!("{:?}", options.initial_centres).into()),
        )
        .with_parameter("maximum_iterations", integer(options.maximum_iterations)?)
        .with_parameter(
            "convergence_tolerance_squared",
            float(options.convergence_tolerance_squared)?,
        );
    Ok(result)
}

#[cfg(test)]
#[path = "kmeans_tests.rs"]
mod tests;
