//! Deterministic Cartesian and circular-feature principal-component analysis.

use molframe_geom::{PeriodicAngle, superpose};
use nalgebra::{DMatrix, SymmetricEigen};
use std::{borrow::Cow, mem::size_of};

use crate::ensemble::{EnsembleGeometryError, generalized_procrustes_mean_source, validate_frames};
use crate::frame_view::{FrameSource, FrameView};
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
    cartesian_pca_source(frames, fit, component_count, memory_limit)
}

/// Runs Cartesian PCA directly over borrowed contiguous frame coordinates.
///
/// Input coordinates are read in place. The `f64` feature matrix is retained
/// as the documented PCA working buffer; frames are never repacked first.
///
/// # Errors
///
/// Returns an error for invalid observations, fitting, controls, or memory requirements.
pub fn cartesian_pca_view(
    frames: FrameView<'_>,
    fit: CartesianFit<'_>,
    component_count: usize,
    memory_limit: usize,
) -> Result<PcaResult, EnsembleGeometryError> {
    cartesian_pca_source(&frames, fit, component_count, memory_limit)
}

fn cartesian_pca_source<S: FrameSource + ?Sized>(
    frames: &S,
    fit: CartesianFit<'_>,
    component_count: usize,
    memory_limit: usize,
) -> Result<PcaResult, EnsembleGeometryError> {
    validate_frames(frames)?;
    let observations = frames.frame_count();
    let features = frames
        .atom_count()
        .checked_mul(3)
        .ok_or(EnsembleGeometryError::MemoryLimit)?;
    let fit_memory = fit_memory(frames.atom_count(), fit)?;
    check_pca_request(
        observations,
        features,
        component_count,
        memory_limit,
        fit_memory,
    )?;
    if let CartesianFit::Reference(reference) = fit
        && reference.len() != frames.atom_count()
    {
        return Err(EnsembleGeometryError::DimensionMismatch);
    }
    let matrix = fitted_matrix(frames, fit, features)?;
    principal_components(matrix, component_count)
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
    let features = first
        .len()
        .checked_mul(2)
        .ok_or(EnsembleGeometryError::MemoryLimit)?;
    check_pca_request(
        observations.len(),
        features,
        component_count,
        memory_limit,
        FitMemory::NONE,
    )?;
    let mut matrix = DMatrix::zeros(observations.len(), features);
    for (row, angles) in observations.iter().enumerate() {
        for (angle, columns) in angles.iter().zip((0..features).step_by(2)) {
            matrix[(row, columns)] = angle.radians().cos();
            matrix[(row, columns + 1)] = angle.radians().sin();
        }
    }
    principal_components(matrix, component_count)
}

fn fitted_matrix<S: FrameSource + ?Sized>(
    frames: &S,
    fit: CartesianFit<'_>,
    features: usize,
) -> Result<DMatrix<f64>, EnsembleGeometryError> {
    let reference = match fit {
        CartesianFit::None => None,
        CartesianFit::Reference(reference) => Some(Cow::Borrowed(reference)),
        CartesianFit::IterativeMean {
            tolerance,
            max_iterations,
        } => Some(Cow::Owned(generalized_procrustes_mean_source(
            frames,
            tolerance,
            max_iterations,
        )?)),
    };
    let mut matrix = DMatrix::zeros(frames.frame_count(), features);
    for index in 0..frames.frame_count() {
        let frame = frames
            .frame(index)
            .ok_or(EnsembleGeometryError::DimensionMismatch)?;
        match reference.as_deref() {
            Some(reference) => {
                let fit = superpose(frame, reference)
                    .map_err(|_| EnsembleGeometryError::DegenerateFit)?;
                for (atom, &point) in frame.iter().enumerate() {
                    write_point(&mut matrix, index, atom, fit.transform.apply(point));
                }
            }
            None => {
                for (atom, &point) in frame.iter().enumerate() {
                    write_point(&mut matrix, index, atom, point);
                }
            }
        }
    }
    Ok(matrix)
}

fn write_point(matrix: &mut DMatrix<f64>, row: usize, atom: usize, point: [f32; 3]) {
    let column = atom * 3;
    matrix[(row, column)] = f64::from(point[0]);
    matrix[(row, column + 1)] = f64::from(point[1]);
    matrix[(row, column + 2)] = f64::from(point[2]);
}

fn principal_components(
    mut centred: DMatrix<f64>,
    component_count: usize,
) -> Result<PcaResult, EnsembleGeometryError> {
    let observations = centred.nrows();
    let features = centred.ncols();
    let mut mean = vec![0.0; features];
    let inverse_count = f64_from_usize(observations)
        .ok_or(EnsembleGeometryError::InvalidParameter)?
        .recip();
    for (column, mean_value) in mean.iter_mut().enumerate() {
        let average = centred.column(column).iter().copied().sum::<f64>() * inverse_count;
        *mean_value = average;
        for value in centred.column_mut(column).iter_mut() {
            *value -= average;
        }
    }
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
    // `tr_mul` forms the product without materialising the transpose, which
    // for a feature-major matrix is a full second copy of the data.
    let covariance = matrix.tr_mul(matrix) / divisor;
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
    let decomposition = SymmetricEigen::new(gram_matrix(matrix, divisor));
    let indices = sorted_indices(decomposition.eigenvalues.as_slice(), limit);
    let mut values = Vec::with_capacity(indices.len());
    let mut vectors = Vec::with_capacity(indices.len());
    for index in indices {
        let value = decomposition.eigenvalues[index].max(0.0);
        if value <= f64::EPSILON {
            continue;
        }
        let snapshot = decomposition.eigenvectors.column(index);
        // Inside the component loop, so materialising the transpose here would
        // allocate a full copy of the matrix once per component.
        let mut vector = matrix.tr_mul(&snapshot);
        vector /= (divisor * value).sqrt();
        let mut feature_vector: Vec<f64> = vector.iter().copied().collect();
        canonicalize_sign(&mut feature_vector);
        values.push(value);
        vectors.push(feature_vector.into_boxed_slice());
    }
    Ok((values, vectors))
}

/// The observation-by-observation Gram matrix, without a transpose copy.
///
/// `matrix * matrix.transpose()` allocates a second copy of the whole matrix,
/// which on the path this is used for — more features than observations, the
/// usual shape for a trajectory — is three doubles per atom per frame.
///
/// Accumulating one column's outer product at a time needs no copy, and reads
/// each column contiguously: `nalgebra` stores column-major, so the row-wise
/// alternative would stride through memory.
fn gram_matrix(matrix: &DMatrix<f64>, divisor: f64) -> DMatrix<f64> {
    let observations = matrix.nrows();
    let mut gram = DMatrix::zeros(observations, observations);

    for feature in 0..matrix.ncols() {
        let column = matrix.column(feature);
        for row in 0..observations {
            let value = column[row];
            if value == 0.0 {
                continue;
            }
            for other in row..observations {
                gram[(row, other)] += value * column[other];
            }
        }
    }

    for row in 0..observations {
        for other in row..observations {
            let value = gram[(row, other)] / divisor;
            gram[(row, other)] = value;
            gram[(other, row)] = value;
        }
    }
    gram
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

fn check_pca_request(
    observations: usize,
    features: usize,
    component_count: usize,
    memory_limit: usize,
    fit_memory: FitMemory,
) -> Result<(), EnsembleGeometryError> {
    if observations < 2 || features == 0 || component_count == 0 {
        return Err(EnsembleGeometryError::InvalidParameter);
    }
    let components = component_count.min(observations - 1).min(features);
    let square = observations.min(features);
    let data = observations
        .checked_mul(features)
        .ok_or(EnsembleGeometryError::MemoryLimit)?;
    let eigensolver = square
        .checked_mul(square)
        .and_then(|entries| entries.checked_mul(2))
        .ok_or(EnsembleGeometryError::MemoryLimit)?;
    let result = features
        .checked_mul(components)
        .and_then(|component_entries| {
            observations
                .checked_mul(components)
                .and_then(|projection_entries| component_entries.checked_add(projection_entries))
        })
        .and_then(|entries| entries.checked_add(features))
        .and_then(|entries| entries.checked_add(components))
        .ok_or(EnsembleGeometryError::MemoryLimit)?;
    let pca_bytes = data
        .checked_add(eigensolver)
        .and_then(|entries| entries.checked_add(result))
        .and_then(|entries| entries.checked_mul(size_of::<f64>()))
        .ok_or(EnsembleGeometryError::MemoryLimit)?;
    let matrix_bytes = data
        .checked_mul(size_of::<f64>())
        .and_then(|bytes| bytes.checked_add(fit_memory.retained_bytes))
        .ok_or(EnsembleGeometryError::MemoryLimit)?;
    let required = pca_bytes
        .max(matrix_bytes)
        .max(fit_memory.construction_peak);
    if required > memory_limit {
        return Err(EnsembleGeometryError::MemoryLimit);
    }
    Ok(())
}

#[derive(Clone, Copy)]
struct FitMemory {
    retained_bytes: usize,
    construction_peak: usize,
}

impl FitMemory {
    const NONE: Self = Self {
        retained_bytes: 0,
        construction_peak: 0,
    };
}

fn fit_memory(atoms: usize, fit: CartesianFit<'_>) -> Result<FitMemory, EnsembleGeometryError> {
    let CartesianFit::IterativeMean {
        tolerance,
        max_iterations,
    } = fit
    else {
        return Ok(FitMemory::NONE);
    };
    if !tolerance.is_finite() || tolerance <= 0.0 || max_iterations == 0 {
        return Err(EnsembleGeometryError::InvalidParameter);
    }
    let retained_bytes = atoms
        .checked_mul(size_of::<[f32; 3]>())
        .ok_or(EnsembleGeometryError::MemoryLimit)?;
    let accumulator_bytes = atoms
        .checked_mul(size_of::<[f64; 3]>())
        .ok_or(EnsembleGeometryError::MemoryLimit)?;
    let construction_peak = retained_bytes
        .checked_mul(2)
        .and_then(|bytes| bytes.checked_add(accumulator_bytes))
        .ok_or(EnsembleGeometryError::MemoryLimit)?;
    Ok(FitMemory {
        retained_bytes,
        construction_peak,
    })
}

#[cfg(test)]
#[path = "pca_tests.rs"]
mod tests;
