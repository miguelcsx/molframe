//! Harmonic and cluster-population similarity between coordinate ensembles.

use nalgebra::{DMatrix, DVector};

use crate::numeric::f64_from_usize;

/// Explicit controls for regularized harmonic ensemble comparison.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicSimilarityOptions {
    /// Positive value added to every covariance diagonal element.
    pub covariance_regularization: f64,
    /// Allocation ceiling for dense covariance matrices.
    pub memory_limit_bytes: usize,
}

/// Bhattacharyya comparison between two regularized Gaussian ensembles.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicSimilarity {
    /// Mean-separation contribution to the distance.
    pub mean_term: f64,
    /// Covariance-volume contribution to the distance.
    pub covariance_term: f64,
    /// Sum of the two Bhattacharyya distance terms.
    pub distance: f64,
    /// `exp(-distance)`, with one denoting identical Gaussian models.
    pub similarity: f64,
}

/// Why ensemble similarity could not be evaluated.
#[derive(Clone, Copy, Debug, thiserror::Error, PartialEq, Eq)]
pub enum EnsembleSimilarityError {
    /// Both ensembles must be non-empty finite rectangular matrices of equal width.
    #[error("ensembles must be non-empty finite matrices with equal feature width")]
    InvalidEnsembles,
    /// Regularization and allocation controls must be finite and positive.
    #[error("regularization and memory limit must be finite and positive")]
    InvalidOptions,
    /// Dense covariance allocation exceeds the caller's ceiling.
    #[error("harmonic comparison exceeds the requested memory limit")]
    MemoryLimit,
    /// The regularized covariance could not be inverted or had invalid determinant.
    #[error("regularized covariance is not positive definite")]
    SingularCovariance,
    /// Cluster labels are empty, out of range, or use no clusters.
    #[error("cluster labels must be non-empty and inside an explicit shared cluster count")]
    InvalidLabels,
}

/// Compares ensembles as regularized multivariate Gaussian distributions.
///
/// # Errors
///
/// Returns [`EnsembleSimilarityError`] for invalid observations, controls,
/// memory requirements, or covariance factorization.
pub fn harmonic_ensemble_similarity(
    first: &[Vec<f64>],
    second: &[Vec<f64>],
    options: HarmonicSimilarityOptions,
) -> Result<HarmonicSimilarity, EnsembleSimilarityError> {
    let dimensions = validate_ensembles(first, second)?;
    validate_options(dimensions, options)?;
    let (first_mean, first_covariance) = gaussian_model(first, options.covariance_regularization)?;
    let (second_mean, second_covariance) =
        gaussian_model(second, options.covariance_regularization)?;
    let average = (&first_covariance + &second_covariance) * 0.5;
    let Some(inverse) = average.clone().try_inverse() else {
        return Err(EnsembleSimilarityError::SingularCovariance);
    };
    let determinant = average.determinant();
    let first_determinant = first_covariance.determinant();
    let second_determinant = second_covariance.determinant();
    if determinant <= 0.0 || first_determinant <= 0.0 || second_determinant <= 0.0 {
        return Err(EnsembleSimilarityError::SingularCovariance);
    }
    let difference = first_mean - second_mean;
    let mean_term = (difference.transpose() * inverse * difference)[(0, 0)] / 8.0;
    let covariance_term =
        0.5 * (determinant / (first_determinant * second_determinant).sqrt()).ln();
    let distance = (mean_term + covariance_term).max(0.0);
    Ok(HarmonicSimilarity {
        mean_term,
        covariance_term,
        distance,
        similarity: (-distance).exp(),
    })
}

/// Compares two ensemble cluster-population distributions using normalized
/// Jensen-Shannon similarity.
///
/// A value of one means identical populations and zero means disjoint support.
/// Cluster identities must already be aligned by the caller.
///
/// # Errors
///
/// Returns [`EnsembleSimilarityError::InvalidLabels`] for absent/out-of-range labels.
pub fn cluster_population_similarity(
    first: &[usize],
    second: &[usize],
    cluster_count: usize,
) -> Result<f64, EnsembleSimilarityError> {
    if first.is_empty()
        || second.is_empty()
        || cluster_count == 0
        || first
            .iter()
            .chain(second)
            .any(|label| *label >= cluster_count)
    {
        return Err(EnsembleSimilarityError::InvalidLabels);
    }
    let first = populations(first, cluster_count)?;
    let second = populations(second, cluster_count)?;
    let midpoint: Vec<f64> = first
        .iter()
        .zip(&second)
        .map(|(left, right)| (left + right) / 2.0)
        .collect();
    let divergence = f64::midpoint(
        relative_entropy(&first, &midpoint),
        relative_entropy(&second, &midpoint),
    );
    Ok(1.0 - divergence / 2.0_f64.ln())
}

fn validate_ensembles(
    first: &[Vec<f64>],
    second: &[Vec<f64>],
) -> Result<usize, EnsembleSimilarityError> {
    let Some(width) = first.first().map(Vec::len) else {
        return Err(EnsembleSimilarityError::InvalidEnsembles);
    };
    if width == 0
        || second.is_empty()
        || first
            .iter()
            .chain(second)
            .any(|row| row.len() != width || row.iter().any(|value| !value.is_finite()))
    {
        return Err(EnsembleSimilarityError::InvalidEnsembles);
    }
    Ok(width)
}

fn validate_options(
    dimensions: usize,
    options: HarmonicSimilarityOptions,
) -> Result<(), EnsembleSimilarityError> {
    if !options.covariance_regularization.is_finite()
        || options.covariance_regularization <= 0.0
        || options.memory_limit_bytes == 0
    {
        return Err(EnsembleSimilarityError::InvalidOptions);
    }
    let required = dimensions
        .checked_mul(dimensions)
        .and_then(|elements| elements.checked_mul(size_of::<f64>()))
        .and_then(|matrix| matrix.checked_mul(4))
        .ok_or(EnsembleSimilarityError::MemoryLimit)?;
    if required > options.memory_limit_bytes {
        return Err(EnsembleSimilarityError::MemoryLimit);
    }
    Ok(())
}

fn gaussian_model(
    observations: &[Vec<f64>],
    regularization: f64,
) -> Result<(DVector<f64>, DMatrix<f64>), EnsembleSimilarityError> {
    let rows =
        f64_from_usize(observations.len()).ok_or(EnsembleSimilarityError::InvalidEnsembles)?;
    let dimensions = observations[0].len();
    let mean = DVector::from_iterator(
        dimensions,
        (0..dimensions)
            .map(|feature| observations.iter().map(|row| row[feature]).sum::<f64>() / rows),
    );
    let mut covariance = DMatrix::zeros(dimensions, dimensions);
    for observation in observations {
        for row in 0..dimensions {
            for column in row..dimensions {
                covariance[(row, column)] +=
                    (observation[row] - mean[row]) * (observation[column] - mean[column]) / rows;
            }
        }
    }
    for row in 0..dimensions {
        for column in 0..row {
            covariance[(row, column)] = covariance[(column, row)];
        }
        covariance[(row, row)] += regularization;
    }
    Ok((mean, covariance))
}

fn populations(
    labels: &[usize],
    cluster_count: usize,
) -> Result<Vec<f64>, EnsembleSimilarityError> {
    let count = f64_from_usize(labels.len()).ok_or(EnsembleSimilarityError::InvalidLabels)?;
    let mut population = vec![0.0; cluster_count];
    for &label in labels {
        population[label] += count.recip();
    }
    Ok(population)
}

fn relative_entropy(distribution: &[f64], reference: &[f64]) -> f64 {
    distribution
        .iter()
        .zip(reference)
        .filter(|(probability, _)| **probability > 0.0)
        .map(|(probability, reference)| probability * (probability / reference).ln())
        .sum()
}

#[cfg(test)]
#[path = "ensemble_similarity_tests.rs"]
mod tests;
