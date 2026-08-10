//! Deterministic Cartesian and circular-feature principal-component analysis.

use nalgebra::{DMatrix, SymmetricEigen};
use pdbiox_geom::{PeriodicAngle, superpose};

use crate::ensemble::{EnsembleGeometryError, generalized_procrustes_mean, validate_frames};
use crate::numeric::f64_from_usize;

type Eigenvectors = (Vec<f64>, Vec<Box<[f64]>>);

/// How Cartesian observations are aligned before PCA.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CartesianFit<'a> {
    /// Preserve the input Cartesian frame.
    None,
    /// Fit each observation to this explicit atom-mapped reference.
    Reference(&'a [[f32; 3]]),
    /// Fit to an iterative generalized Procrustes mean.
    IterativeMean {
        /// RMSD convergence threshold.
        tolerance: f64,
        /// Maximum number of mean iterations.
        max_iterations: usize,
    },
}

/// PCA eigenvalues, feature-space axes and observation projections.
#[derive(Clone, Debug, PartialEq)]
pub struct PcaResult {
    /// Mean of each fitted feature.
    pub mean: Box<[f64]>,
    /// Eigenvalues in descending order.
    pub eigenvalues: Box<[f64]>,
    /// Feature-space axes, one canonical-sign vector per component.
    pub components: Box<[Box<[f64]>]>,
    /// Observation coordinates, one vector per observation.
    pub projections: Box<[Box<[f64]>]>,
}

/// Runs Cartesian PCA with an explicit fitting policy and allocation ceiling.
///
/// # Errors
///
/// Returns an error for invalid observations, fitting, controls, or memory requirements.
pub fn cartesian_pca(
    frames: &[Vec<[f32; 3]>],
    fit: CartesianFit<'_>,
    component_count: usize,
    memory_limit: usize,
) -> Result<PcaResult, EnsembleGeometryError> {
    validate_frames(frames)?;
    let fitted = fit_frames(frames, fit)?;
    let rows: Vec<Vec<f64>> = fitted
        .iter()
        .map(|frame| {
            frame
                .iter()
                .flat_map(|point| point.iter().map(|value| f64::from(*value)))
                .collect()
        })
        .collect();
    principal_components(&rows, component_count, memory_limit)
}

/// Runs dihedral PCA on `(cos θ, sin θ)` features.
///
/// This embedding is invariant to adding any integer multiple of `2π`.
///
/// # Errors
///
/// Returns an error for inconsistent observations, controls, or memory requirements.
pub fn dihedral_pca(
    observations: &[Vec<PeriodicAngle>],
    component_count: usize,
    memory_limit: usize,
) -> Result<PcaResult, EnsembleGeometryError> {
    let Some(first) = observations.first() else {
        return Err(EnsembleGeometryError::Empty);
    };
    if first.is_empty() || observations.iter().any(|row| row.len() != first.len()) {
        return Err(EnsembleGeometryError::DimensionMismatch);
    }
    let rows: Vec<Vec<f64>> = observations
        .iter()
        .map(|angles| {
            angles
                .iter()
                .flat_map(|angle| [angle.radians().cos(), angle.radians().sin()])
                .collect()
        })
        .collect();
    principal_components(&rows, component_count, memory_limit)
}

fn fit_frames(
    frames: &[Vec<[f32; 3]>],
    fit: CartesianFit<'_>,
) -> Result<Vec<Vec<[f32; 3]>>, EnsembleGeometryError> {
    let reference = match fit {
        CartesianFit::None => return Ok(frames.to_vec()),
        CartesianFit::Reference(reference) => {
            if reference.len() != frames[0].len() {
                return Err(EnsembleGeometryError::DimensionMismatch);
            }
            reference.to_vec()
        }
        CartesianFit::IterativeMean {
            tolerance,
            max_iterations,
        } => generalized_procrustes_mean(frames, tolerance, max_iterations)?,
    };
    frames
        .iter()
        .map(|frame| {
            let fit =
                superpose(frame, &reference).map_err(|_| EnsembleGeometryError::DegenerateFit)?;
            Ok(frame
                .iter()
                .map(|&point| fit.transform.apply(point))
                .collect())
        })
        .collect()
}

fn principal_components(
    rows: &[Vec<f64>],
    component_count: usize,
    memory_limit: usize,
) -> Result<PcaResult, EnsembleGeometryError> {
    if rows.len() < 2 || component_count == 0 {
        return Err(EnsembleGeometryError::InvalidParameter);
    }
    let features = rows[0].len();
    if features == 0 || rows.iter().any(|row| row.len() != features) {
        return Err(EnsembleGeometryError::DimensionMismatch);
    }
    let observations = rows.len();
    check_pca_memory(observations, features, memory_limit)?;
    let mut mean = vec![0.0; features];
    for row in rows {
        for (average, value) in mean.iter_mut().zip(row) {
            *average += value;
        }
    }
    let inverse_count = f64_from_usize(observations)
        .ok_or(EnsembleGeometryError::InvalidParameter)?
        .recip();
    for average in &mut mean {
        *average *= inverse_count;
    }
    let centred_values: Vec<f64> = rows
        .iter()
        .flat_map(|row| {
            row.iter()
                .zip(&mean)
                .map(|(value, average)| value - average)
        })
        .collect();
    let centred = DMatrix::from_row_slice(observations, features, &centred_values);
    let limit = component_count.min(observations - 1).min(features);
    let (eigenvalues, components) = if features <= observations {
        feature_eigenvectors(&centred, limit)?
    } else {
        snapshot_eigenvectors(&centred, limit)?
    };
    let projections = project(&centred, &components);
    Ok(PcaResult {
        mean: mean.into_boxed_slice(),
        eigenvalues: eigenvalues.into_boxed_slice(),
        components: components.into_boxed_slice(),
        projections,
    })
}

fn feature_eigenvectors(
    matrix: &DMatrix<f64>,
    limit: usize,
) -> Result<Eigenvectors, EnsembleGeometryError> {
    let divisor = observation_divisor(matrix)?;
    let covariance = matrix.transpose() * matrix / divisor;
    Ok(selected_eigenvectors(
        &SymmetricEigen::new(covariance),
        limit,
    ))
}

fn snapshot_eigenvectors(
    matrix: &DMatrix<f64>,
    limit: usize,
) -> Result<Eigenvectors, EnsembleGeometryError> {
    let divisor = observation_divisor(matrix)?;
    let gram = matrix * matrix.transpose() / divisor;
    let decomposition = SymmetricEigen::new(gram);
    let indices = sorted_indices(decomposition.eigenvalues.as_slice(), limit);
    let mut values = Vec::with_capacity(indices.len());
    let mut vectors = Vec::with_capacity(indices.len());
    for index in indices {
        let value = decomposition.eigenvalues[index].max(0.0);
        if value <= f64::EPSILON {
            continue;
        }
        let snapshot = decomposition.eigenvectors.column(index);
        let mut vector = matrix.transpose() * snapshot;
        vector /= (divisor * value).sqrt();
        let mut feature_vector: Vec<f64> = vector.iter().copied().collect();
        canonicalize_sign(&mut feature_vector);
        values.push(value);
        vectors.push(feature_vector.into_boxed_slice());
    }
    Ok((values, vectors))
}

fn observation_divisor(matrix: &DMatrix<f64>) -> Result<f64, EnsembleGeometryError> {
    matrix
        .nrows()
        .checked_sub(1)
        .and_then(f64_from_usize)
        .ok_or(EnsembleGeometryError::InvalidParameter)
}

fn selected_eigenvectors(
    decomposition: &SymmetricEigen<f64, nalgebra::Dyn>,
    limit: usize,
) -> (Vec<f64>, Vec<Box<[f64]>>) {
    let indices = sorted_indices(decomposition.eigenvalues.as_slice(), limit);
    let mut values = Vec::with_capacity(indices.len());
    let mut vectors = Vec::with_capacity(indices.len());
    for index in indices {
        let value = decomposition.eigenvalues[index].max(0.0);
        let mut vector: Vec<f64> = decomposition
            .eigenvectors
            .column(index)
            .iter()
            .copied()
            .collect();
        canonicalize_sign(&mut vector);
        values.push(value);
        vectors.push(vector.into_boxed_slice());
    }
    (values, vectors)
}

fn sorted_indices(values: &[f64], limit: usize) -> Vec<usize> {
    let mut indices: Vec<usize> = (0..values.len()).collect();
    indices.sort_by(|&left, &right| {
        values[right]
            .total_cmp(&values[left])
            .then(left.cmp(&right))
    });
    indices.truncate(limit);
    indices
}

fn canonicalize_sign(vector: &mut [f64]) {
    let pivot = vector
        .iter()
        .enumerate()
        .max_by(|(left_index, left), (right_index, right)| {
            left.abs()
                .total_cmp(&right.abs())
                .then(right_index.cmp(left_index))
        });
    if pivot.is_some_and(|(_, value)| *value < 0.0) {
        for value in vector {
            *value = -*value;
        }
    }
}

fn project(matrix: &DMatrix<f64>, components: &[Box<[f64]>]) -> Box<[Box<[f64]>]> {
    (0..matrix.nrows())
        .map(|row| {
            components
                .iter()
                .map(|component| {
                    component
                        .iter()
                        .enumerate()
                        .map(|(column, weight)| matrix[(row, column)] * weight)
                        .sum()
                })
                .collect::<Vec<f64>>()
                .into_boxed_slice()
        })
        .collect()
}

fn check_pca_memory(
    observations: usize,
    features: usize,
    limit: usize,
) -> Result<(), EnsembleGeometryError> {
    let square = observations.min(features);
    let entries = observations
        .checked_mul(features)
        .and_then(|data| {
            square
                .checked_mul(square)
                .and_then(|covariance| data.checked_add(covariance))
        })
        .and_then(|count| count.checked_mul(size_of::<f64>()))
        .ok_or(EnsembleGeometryError::MemoryLimit)?;
    if entries > limit {
        return Err(EnsembleGeometryError::MemoryLimit);
    }
    Ok(())
}

#[cfg(test)]
#[path = "pca_tests.rs"]
mod tests;
