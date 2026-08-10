//! Deterministic k-means with caller-selected initial centres.

use std::collections::BTreeSet;

use crate::numeric::f64_from_usize;

/// Explicit k-means convergence and initialization policy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KMeansOptions<'a> {
    /// Observation indices copied as initial centres, one per cluster.
    pub initial_centres: &'a [usize],
    /// Maximum assignment/update iterations.
    pub maximum_iterations: usize,
    /// Maximum squared centre displacement accepted as converged.
    pub convergence_tolerance_squared: f64,
}

/// A converged k-means partition.
#[derive(Clone, Debug, PartialEq)]
pub struct KMeans {
    /// Cluster index per observation.
    pub labels: Vec<usize>,
    /// Arithmetic mean feature vector per cluster.
    pub centres: Vec<Vec<f64>>,
    /// Sum of squared distances to assigned centres.
    pub inertia: f64,
    /// Iterations required to converge.
    pub iterations: usize,
}

/// Why k-means could not be evaluated.
#[derive(Clone, Copy, Debug, thiserror::Error, PartialEq, Eq)]
pub enum KMeansError {
    /// Observations must be a non-empty finite rectangular matrix.
    #[error("observations must be a non-empty finite rectangular matrix")]
    InvalidObservations,
    /// Initialization and convergence controls are invalid.
    #[error(
        "initial centres must be unique and in range; iteration and tolerance controls must be valid"
    )]
    InvalidOptions,
    /// A cluster lost all observations during an update.
    #[error("cluster {0} became empty")]
    EmptyCluster(usize),
    /// The explicit iteration ceiling was reached before convergence.
    #[error("k-means did not converge within the requested iteration ceiling")]
    DidNotConverge,
}

/// Clusters feature vectors using explicit deterministic initialization.
///
/// No random seed, implicit `k`, or k-means++ heuristic is hidden in this API.
/// Equal-distance assignments choose the lower cluster index.
///
/// # Errors
///
/// Returns [`KMeansError`] for invalid data/options, an empty cluster, or non-convergence.
pub fn kmeans(
    observations: &[Vec<f64>],
    options: KMeansOptions<'_>,
) -> Result<KMeans, KMeansError> {
    let dimensions = validate(observations, options)?;
    let mut centres: Vec<Vec<f64>> = options
        .initial_centres
        .iter()
        .map(|index| observations[*index].clone())
        .collect();
    let mut labels = vec![0; observations.len()];
    for iteration in 1..=options.maximum_iterations {
        assign(observations, &centres, &mut labels);
        let next = update(observations, &labels, centres.len(), dimensions)?;
        let movement = centres
            .iter()
            .zip(&next)
            .map(|(left, right)| squared_distance(left, right))
            .fold(0.0, f64::max);
        centres = next;
        if movement <= options.convergence_tolerance_squared {
            let inertia = observations
                .iter()
                .zip(&labels)
                .map(|(observation, label)| squared_distance(observation, &centres[*label]))
                .sum();
            return Ok(KMeans {
                labels,
                centres,
                inertia,
                iterations: iteration,
            });
        }
    }
    Err(KMeansError::DidNotConverge)
}

fn validate(observations: &[Vec<f64>], options: KMeansOptions<'_>) -> Result<usize, KMeansError> {
    let Some(first) = observations.first() else {
        return Err(KMeansError::InvalidObservations);
    };
    if first.is_empty()
        || observations.iter().any(|row| row.len() != first.len())
        || observations
            .iter()
            .flatten()
            .any(|value| !value.is_finite())
    {
        return Err(KMeansError::InvalidObservations);
    }
    let unique: BTreeSet<usize> = options.initial_centres.iter().copied().collect();
    if options.initial_centres.is_empty()
        || unique.len() != options.initial_centres.len()
        || unique.iter().any(|index| *index >= observations.len())
        || options.maximum_iterations == 0
        || !options.convergence_tolerance_squared.is_finite()
        || options.convergence_tolerance_squared < 0.0
    {
        return Err(KMeansError::InvalidOptions);
    }
    Ok(first.len())
}

fn assign(observations: &[Vec<f64>], centres: &[Vec<f64>], labels: &mut [usize]) {
    for (observation, label) in observations.iter().zip(labels) {
        *label = centres
            .iter()
            .enumerate()
            .map(|(cluster, centre)| (squared_distance(observation, centre), cluster))
            .min_by(|left, right| left.0.total_cmp(&right.0).then(left.1.cmp(&right.1)))
            .map_or(0, |(_, cluster)| cluster);
    }
}

fn update(
    observations: &[Vec<f64>],
    labels: &[usize],
    cluster_count: usize,
    dimensions: usize,
) -> Result<Vec<Vec<f64>>, KMeansError> {
    let mut centres = vec![vec![0.0; dimensions]; cluster_count];
    let mut counts = vec![0_usize; cluster_count];
    for (observation, &label) in observations.iter().zip(labels) {
        counts[label] += 1;
        for (sum, value) in centres[label].iter_mut().zip(observation) {
            *sum += value;
        }
    }
    for (cluster, (centre, count)) in centres.iter_mut().zip(counts).enumerate() {
        if count == 0 {
            return Err(KMeansError::EmptyCluster(cluster));
        }
        let divisor = f64_from_usize(count).ok_or(KMeansError::InvalidObservations)?;
        for value in centre {
            *value /= divisor;
        }
    }
    Ok(centres)
}

fn squared_distance(left: &[f64], right: &[f64]) -> f64 {
    left.iter()
        .zip(right)
        .map(|(left, right)| (left - right).powi(2))
        .sum()
}

#[cfg(test)]
#[path = "kmeans_tests.rs"]
mod tests;
