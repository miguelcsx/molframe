//! Pairwise geometry and generalized Procrustes summaries for ensembles.

use pdbiox_geom::{PeriodicAngle, TorusMetric, rmsd, superpose};
use thiserror::Error;

use crate::frame_view::{FrameSource, FrameView};
use crate::numeric::{f32_triplet, f64_from_usize};

/// Default ceiling for a public pairwise allocation.
pub const DEFAULT_PAIRWISE_MEMORY_LIMIT: usize = 512 * 1024 * 1024;
/// Default convergence tolerance for a generalized Procrustes mean.
pub const DEFAULT_PROCRUSTES_TOLERANCE: f64 = 1.0e-6;
/// Default iteration ceiling for a generalized Procrustes mean.
pub const DEFAULT_PROCRUSTES_MAXIMUM_ITERATIONS: usize = 100;

/// Why an ensemble calculation could not be completed.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Error)]
pub enum EnsembleGeometryError {
    /// No observations were supplied.
    #[error("the ensemble is empty")]
    Empty,
    /// Observations do not share one feature or atom count.
    #[error("ensemble observations have different dimensions")]
    DimensionMismatch,
    /// A rigid fit is undefined for the supplied coordinates.
    #[error("a rigid fit is undefined")]
    DegenerateFit,
    /// A requested parameter is not finite or positive.
    #[error("a parameter must be finite and positive")]
    InvalidParameter,
    /// The requested dense allocation exceeds its explicit ceiling.
    #[error("the calculation exceeds its memory limit")]
    MemoryLimit,
    /// Iteration stopped before satisfying the convergence tolerance.
    #[error("the iterative mean did not converge")]
    DidNotConverge,
}

/// Dense symmetric pairwise distances in row-major order.
#[derive(Clone, Debug, PartialEq)]
pub struct EnsembleDistanceMatrix {
    /// Number of observations on each matrix axis.
    pub size: usize,
    /// Row-major distances, including a zero diagonal.
    pub values: Box<[f64]>,
}

/// Alignment policy for RMSD against one trajectory frame.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FrameAlignment {
    /// Measure coordinates in their current frame.
    None,
    /// Rigidly fit each frame onto the reference before measuring.
    Rigid,
}

/// Computes one RMSD per frame against a selected reference frame.
///
/// # Errors
///
/// Returns an error for an empty trajectory, an absent reference, changing atom
/// counts, or a frame that cannot define the requested rigid fit.
pub fn rmsd_to_reference(
    frames: &[crate::Timestep],
    reference: usize,
    alignment: FrameAlignment,
) -> Result<Vec<f64>, EnsembleGeometryError> {
    rmsd_to_reference_source(frames, reference, alignment)
}

/// Computes one RMSD per borrowed frame against a selected reference frame.
///
/// The view is consumed in place without materialising owned trajectory frames.
///
/// # Errors
///
/// Returns an error for empty, malformed, incompatible, or degenerate frames.
pub fn rmsd_to_reference_view(
    frames: FrameView<'_>,
    reference: usize,
    alignment: FrameAlignment,
) -> Result<Vec<f64>, EnsembleGeometryError> {
    rmsd_to_reference_source(&frames, reference, alignment)
}

fn rmsd_to_reference_source<S: FrameSource + ?Sized>(
    frames: &S,
    reference: usize,
    alignment: FrameAlignment,
) -> Result<Vec<f64>, EnsembleGeometryError> {
    if frames.frame_count() == 0 {
        return Err(EnsembleGeometryError::Empty);
    }
    let target = frames
        .frame(reference)
        .ok_or(EnsembleGeometryError::InvalidParameter)?;
    let mut values = Vec::with_capacity(frames.frame_count());
    for index in 0..frames.frame_count() {
        let frame = frames
            .frame(index)
            .ok_or(EnsembleGeometryError::DimensionMismatch)?;
        if frame.len() != target.len() {
            return Err(EnsembleGeometryError::DimensionMismatch);
        }
        let value = match alignment {
            FrameAlignment::None => {
                rmsd(frame, target).map_err(|_| EnsembleGeometryError::DimensionMismatch)?
            }
            FrameAlignment::Rigid => superpose(frame, target)
                .map(|fit| fit.rmsd)
                .map_err(|_| EnsembleGeometryError::DegenerateFit)?,
        };
        values.push(value);
    }
    Ok(values)
}

impl EnsembleDistanceMatrix {
    /// Returns a matrix entry.
    #[must_use]
    pub fn get(&self, row: usize, column: usize) -> Option<f64> {
        row.checked_mul(self.size)
            .and_then(|offset| offset.checked_add(column))
            .and_then(|index| self.values.get(index).copied())
    }
}

/// Computes all fitted RMSD pairs after validating the allocation ceiling.
///
/// # Errors
///
/// Returns an error for inconsistent or degenerate frames, or an oversized allocation.
pub fn pairwise_fitted_rmsd(
    frames: &[Vec<[f32; 3]>],
    memory_limit: usize,
) -> Result<EnsembleDistanceMatrix, EnsembleGeometryError> {
    pairwise_fitted_rmsd_source(frames, memory_limit)
}

/// Computes all fitted RMSD pairs directly from borrowed contiguous frames.
///
/// The result allocation is bounded by `memory_limit`; input coordinates remain
/// borrowed for the full computation.
///
/// # Errors
///
/// Returns an error for invalid frames, an oversized result, or a degenerate fit.
pub fn pairwise_fitted_rmsd_view(
    frames: FrameView<'_>,
    memory_limit: usize,
) -> Result<EnsembleDistanceMatrix, EnsembleGeometryError> {
    pairwise_fitted_rmsd_source(&frames, memory_limit)
}

fn pairwise_fitted_rmsd_source<S: FrameSource + ?Sized>(
    frames: &S,
    memory_limit: usize,
) -> Result<EnsembleDistanceMatrix, EnsembleGeometryError> {
    validate_frames(frames)?;
    let frame_count = frames.frame_count();
    let mut matrix = allocate_square(frame_count, memory_limit)?;
    for left in 0..frame_count {
        let left_frame = frames
            .frame(left)
            .ok_or(EnsembleGeometryError::DimensionMismatch)?;
        for right in (left + 1)..frame_count {
            let right_frame = frames
                .frame(right)
                .ok_or(EnsembleGeometryError::DimensionMismatch)?;
            let fit = superpose(left_frame, right_frame)
                .map_err(|_| EnsembleGeometryError::DegenerateFit)?;
            set_pair(&mut matrix, frame_count, left, right, fit.rmsd);
        }
    }
    Ok(EnsembleDistanceMatrix {
        size: frame_count,
        values: matrix.into_boxed_slice(),
    })
}

/// Computes all torus-distance pairs after validating dimensions and memory.
///
/// # Errors
///
/// Returns an error for incompatible dimensions or an oversized allocation.
pub fn pairwise_torus_distance(
    points: &[Vec<PeriodicAngle>],
    metric: &TorusMetric,
    memory_limit: usize,
) -> Result<EnsembleDistanceMatrix, EnsembleGeometryError> {
    if points.is_empty() {
        return Err(EnsembleGeometryError::Empty);
    }
    let mut matrix = allocate_square(points.len(), memory_limit)?;
    for left in 0..points.len() {
        for right in (left + 1)..points.len() {
            let distance = metric
                .distance(&points[left], &points[right])
                .map_err(|_| EnsembleGeometryError::DimensionMismatch)?;
            set_pair(&mut matrix, points.len(), left, right, distance);
        }
    }
    Ok(EnsembleDistanceMatrix {
        size: points.len(),
        values: matrix.into_boxed_slice(),
    })
}

/// Finds an iterative generalized Procrustes mean under explicit convergence controls.
///
/// # Errors
///
/// Returns an error for invalid frames, controls, fits, or failed convergence.
pub fn generalized_procrustes_mean(
    frames: &[Vec<[f32; 3]>],
    tolerance: f64,
    max_iterations: usize,
) -> Result<Vec<[f32; 3]>, EnsembleGeometryError> {
    generalized_procrustes_mean_source(frames, tolerance, max_iterations)
}

/// Finds a generalized Procrustes mean directly from borrowed contiguous frames.
///
/// The evolving mean and per-atom accumulator are owned working buffers; input
/// coordinates are never repacked into per-frame vectors.
///
/// # Errors
///
/// Returns an error for invalid frames, controls, fits, or failed convergence.
pub fn generalized_procrustes_mean_view(
    frames: FrameView<'_>,
    tolerance: f64,
    max_iterations: usize,
) -> Result<Vec<[f32; 3]>, EnsembleGeometryError> {
    generalized_procrustes_mean_source(&frames, tolerance, max_iterations)
}

pub(crate) fn generalized_procrustes_mean_source<S: FrameSource + ?Sized>(
    frames: &S,
    tolerance: f64,
    max_iterations: usize,
) -> Result<Vec<[f32; 3]>, EnsembleGeometryError> {
    validate_frames(frames)?;
    if !tolerance.is_finite() || tolerance <= 0.0 || max_iterations == 0 {
        return Err(EnsembleGeometryError::InvalidParameter);
    }
    let first = frames
        .frame(0)
        .ok_or(EnsembleGeometryError::DimensionMismatch)?;
    let mut mean = first.to_vec();
    for _ in 0..max_iterations {
        let mut sum = vec![[0.0_f64; 3]; mean.len()];
        for index in 0..frames.frame_count() {
            let frame = frames
                .frame(index)
                .ok_or(EnsembleGeometryError::DimensionMismatch)?;
            let fit = superpose(frame, &mean).map_err(|_| EnsembleGeometryError::DegenerateFit)?;
            for (accumulator, &point) in sum.iter_mut().zip(frame) {
                let transformed = fit.transform.apply(point);
                for axis in 0..3 {
                    accumulator[axis] += f64::from(transformed[axis]);
                }
            }
        }
        let inverse = f64_from_usize(frames.frame_count())
            .ok_or(EnsembleGeometryError::InvalidParameter)?
            .recip();
        let next: Vec<_> = sum
            .into_iter()
            .map(|point| {
                f32_triplet(point.map(|value| value * inverse))
                    .ok_or(EnsembleGeometryError::InvalidParameter)
            })
            .collect::<Result<_, _>>()?;
        let change = rmsd(&mean, &next).map_err(|_| EnsembleGeometryError::DegenerateFit)?;
        mean = next;
        if change <= tolerance {
            return Ok(mean);
        }
    }
    Err(EnsembleGeometryError::DidNotConverge)
}

pub(crate) fn validate_frames<S: FrameSource + ?Sized>(
    frames: &S,
) -> Result<usize, EnsembleGeometryError> {
    if frames.frame_count() == 0 {
        return Err(EnsembleGeometryError::Empty);
    }
    let atom_count = frames.atom_count();
    if atom_count < 3 {
        return Err(EnsembleGeometryError::DimensionMismatch);
    }
    for index in 0..frames.frame_count() {
        let Some(frame) = frames.frame(index) else {
            return Err(EnsembleGeometryError::DimensionMismatch);
        };
        if frame.len() != atom_count {
            return Err(EnsembleGeometryError::DimensionMismatch);
        }
    }
    Ok(atom_count)
}

pub(crate) fn allocate_square(
    size: usize,
    memory_limit: usize,
) -> Result<Vec<f64>, EnsembleGeometryError> {
    let bytes = size
        .checked_mul(size)
        .and_then(|count| count.checked_mul(size_of::<f64>()))
        .ok_or(EnsembleGeometryError::MemoryLimit)?;
    if bytes > memory_limit {
        return Err(EnsembleGeometryError::MemoryLimit);
    }
    Ok(vec![0.0; size * size])
}

fn set_pair(matrix: &mut [f64], size: usize, left: usize, right: usize, value: f64) {
    matrix[left * size + right] = value;
    matrix[right * size + left] = value;
}

#[cfg(test)]
#[path = "ensemble_tests.rs"]
mod tests;
