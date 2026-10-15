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

/// Borrowed row-major `f64` observations for zero-copy numerical callers.
///
/// This is kept private to the implementation; [`kmeans_view`] is the public
/// shape-checked entry point so FFI callers do not need to materialise nested
/// `Vec`s merely to satisfy the owned Rust convenience API.
#[derive(Clone, Copy, Debug)]
struct FlatObservations<'a> {
    values: &'a [f64],
    rows: usize,
    columns: usize,
}

impl<'a> FlatObservations<'a> {
    fn new(values: &'a [f64], rows: usize, columns: usize) -> Option<Self> {
        let expected = rows.checked_mul(columns)?;
        (expected == values.len()).then_some(Self {
            values,
            rows,
            columns,
        })
    }
}

trait ObservationSource {
    fn row_count(&self) -> usize;
    fn column_count(&self) -> usize;
    fn row(&self, index: usize) -> Option<&[f64]>;
}

impl ObservationSource for [Vec<f64>] {
    fn row_count(&self) -> usize {
        self.len()
    }

    fn column_count(&self) -> usize {
        self.first().map_or(0, Vec::len)
    }

    fn row(&self, index: usize) -> Option<&[f64]> {
        self.get(index).map(Vec::as_slice)
    }
}

impl ObservationSource for FlatObservations<'_> {
    fn row_count(&self) -> usize {
        self.rows
    }

    fn column_count(&self) -> usize {
        self.columns
    }

    fn row(&self, index: usize) -> Option<&[f64]> {
        let start = index.checked_mul(self.columns)?;
        let end = start.checked_add(self.columns)?;
        self.values.get(start..end)
    }
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
    kmeans_source(observations, options)
}

/// Clusters C-contiguous row-major observations without repacking them.
///
/// `values` must contain exactly `rows * columns` finite values in row-major
/// order.  The algorithm only allocates its labels, centres, and result; input
/// observations remain borrowed for the full kernel execution.
///
/// # Errors
///
/// Returns [`KMeansError::InvalidObservations`] when dimensions overflow, the
/// shape does not match the slice, or data are not a finite non-empty matrix.
pub fn kmeans_view(
    values: &[f64],
    rows: usize,
    columns: usize,
    options: KMeansOptions<'_>,
) -> Result<KMeans, KMeansError> {
    let observations =
        FlatObservations::new(values, rows, columns).ok_or(KMeansError::InvalidObservations)?;
    kmeans_source(&observations, options)
}

fn kmeans_source<S: ObservationSource + ?Sized>(
    observations: &S,
    options: KMeansOptions<'_>,
) -> Result<KMeans, KMeansError> {
    let dimensions = validate(observations, options)?;
    let mut centres: Vec<Vec<f64>> = options
        .initial_centres
        .iter()
        .map(|index| {
            observations
                .row(*index)
                .map(<[f64]>::to_vec)
                .ok_or(KMeansError::InvalidObservations)
        })
        .collect::<Result<_, _>>()?;
    let mut labels = vec![0; observations.row_count()];
    for iteration in 1..=options.maximum_iterations {
        assign(observations, &centres, &mut labels)?;
        let next = update(observations, &labels, centres.len(), dimensions)?;
        let movement = centres
            .iter()
            .zip(&next)
            .map(|(left, right)| squared_distance(left, right))
            .fold(0.0, f64::max);
        centres = next;
        if movement <= options.convergence_tolerance_squared {
            let inertia = (0..observations.row_count()).try_fold(0.0, |sum, index| {
                let observation = observations
                    .row(index)
                    .ok_or(KMeansError::InvalidObservations)?;
                Ok(sum + squared_distance(observation, &centres[labels[index]]))
            })?;
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

fn validate<S: ObservationSource + ?Sized>(
    observations: &S,
    options: KMeansOptions<'_>,
) -> Result<usize, KMeansError> {
    let dimensions = observations.column_count();
    if observations.row_count() == 0 || dimensions == 0 {
        return Err(KMeansError::InvalidObservations);
    }
    for index in 0..observations.row_count() {
        let row = observations
            .row(index)
            .ok_or(KMeansError::InvalidObservations)?;
        if row.len() != dimensions || row.iter().any(|value| !value.is_finite()) {
            return Err(KMeansError::InvalidObservations);
        }
    }
    let unique: BTreeSet<usize> = options.initial_centres.iter().copied().collect();
    if options.initial_centres.is_empty()
        || unique.len() != options.initial_centres.len()
        || unique
            .iter()
            .any(|index| *index >= observations.row_count())
        || options.maximum_iterations == 0
        || !options.convergence_tolerance_squared.is_finite()
        || options.convergence_tolerance_squared < 0.0
    {
        return Err(KMeansError::InvalidOptions);
    }
    Ok(dimensions)
}

fn assign<S: ObservationSource + ?Sized>(
    observations: &S,
    centres: &[Vec<f64>],
    labels: &mut [usize],
) -> Result<(), KMeansError> {
    if labels.len() != observations.row_count() {
        return Err(KMeansError::InvalidObservations);
    }
    for (index, label) in labels.iter_mut().enumerate() {
        let observation = observations
            .row(index)
            .ok_or(KMeansError::InvalidObservations)?;
        *label = centres
            .iter()
            .enumerate()
            .map(|(cluster, centre)| (squared_distance(observation, centre), cluster))
            .min_by(|left, right| left.0.total_cmp(&right.0).then(left.1.cmp(&right.1)))
            .map_or(0, |(_, cluster)| cluster);
    }
    Ok(())
}

fn update<S: ObservationSource + ?Sized>(
    observations: &S,
    labels: &[usize],
    cluster_count: usize,
    dimensions: usize,
) -> Result<Vec<Vec<f64>>, KMeansError> {
    let mut centres = vec![vec![0.0; dimensions]; cluster_count];
    let mut counts = vec![0_usize; cluster_count];
    for (index, &label) in labels.iter().enumerate() {
        let observation = observations
            .row(index)
            .ok_or(KMeansError::InvalidObservations)?;
        if label >= cluster_count {
            return Err(KMeansError::InvalidObservations);
        }
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
