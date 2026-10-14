//! Deterministic Cartesian and circular-feature principal-component analysis.
//!
//! For `n` observations of `p = 3N` features the observation matrix is held
//! once, row-major, in `f64` (`8 n p` bytes) and is never centred in place or
//! transposed; products centre one block at a time on the worker pool.
//!
//! | route | chosen when | time | extra memory |
//! |---|---|---|---|
//! | covariance | `p <= n` | `O(n p^2 + p^3)` | `p^2` per in-flight block, `2 p^2` solver |
//! | snapshot (Gram) | `p > n` | `O(n^2 p + n^3)` | `n^2` per in-flight block, `2 n^2` solver |
//!
//! The snapshot route recovers feature-space axes as `X^T v / sqrt((n-1) l)`
//! in one parallel product, so the `p x p` covariance never exists when
//! `n < p`. Eigenvalues agree with a dense covariance reference to a relative
//! `1e-9` and axes to `1e-7` (unit vectors, after sign canonicalisation) on
//! well-separated spectra; see `pca_reference_tests.rs`.

use crate::ensemble::{EnsembleGeometryError, validate_frames};
use crate::frame_view::{FrameSource, FrameView};
use crate::numeric::f64_from_usize;
use crate::pca_fit::{column_means, fitted_observations};
use crate::pca_linalg::{Observations, back_project, covariance, gram, project, top_eigenpairs};
use crate::pca_memory::{FitMemory, check_pca_request, fit_memory};
use faer::Mat;
use molframe_core::ExecutionContext;
use molframe_geom::PeriodicAngle;

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

fn cartesian_pca_source<S: FrameSource + Sync + ?Sized>(
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
    let context = ExecutionContext::default();
    let fitted = fitted_observations(frames, fit, features, &context)?;
    principal_components(
        &fitted.data,
        &fitted.mean,
        observations,
        features,
        component_count,
        &context,
    )
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
    let mut data = Vec::with_capacity(observations.len() * features);
    for angles in observations {
        for angle in angles {
            data.push(angle.radians().cos());
            data.push(angle.radians().sin());
        }
    }
    let mean = column_means(&data, features)?;
    principal_components(
        &data,
        &mean,
        observations.len(),
        features,
        component_count,
        &ExecutionContext::default(),
    )
}

fn principal_components(
    data: &[f64],
    mean: &[f64],
    observations: usize,
    features: usize,
    component_count: usize,
    context: &ExecutionContext,
) -> Result<PcaResult, EnsembleGeometryError> {
    let matrix = Observations::new(data, observations, features, mean)?;
    let limit = component_count.min(observations - 1).min(features);
    let divisor =
        f64_from_usize(observations - 1).ok_or(EnsembleGeometryError::InvalidParameter)?;
    let (eigenvalues, mut components) = if features <= observations {
        feature_eigenvectors(matrix, divisor, limit, context)?
    } else {
        snapshot_eigenvectors(matrix, divisor, limit, context)?
    };
    canonicalize_signs(&mut components);
    let projections = project(matrix, components.as_ref(), context)?;
    Ok(PcaResult {
        mean: mean.into(),
        eigenvalues: eigenvalues.into_boxed_slice(),
        components: (0..components.ncols())
            .map(|column| (0..features).map(|row| components[(row, column)]).collect())
            .collect(),
        projections,
    })
}

fn feature_eigenvectors(
    matrix: Observations<'_>,
    divisor: f64,
    limit: usize,
    context: &ExecutionContext,
) -> Result<(Vec<f64>, Mat<f64>), EnsembleGeometryError> {
    let (values, vectors) = top_eigenpairs(&covariance(matrix, divisor, context)?, limit)?;
    Ok((
        values.into_iter().map(|value| value.max(0.0)).collect(),
        vectors,
    ))
}

fn snapshot_eigenvectors(
    matrix: Observations<'_>,
    divisor: f64,
    limit: usize,
    context: &ExecutionContext,
) -> Result<(Vec<f64>, Mat<f64>), EnsembleGeometryError> {
    let (values, snapshots) = top_eigenpairs(&gram(matrix, divisor, context)?, limit)?;
    let values: Vec<f64> = values
        .into_iter()
        .map(|value| value.max(0.0))
        .take_while(|&value| value > f64::EPSILON)
        .collect();
    let kept = snapshots.as_ref().subcols(0, values.len());
    let mut axes = back_project(matrix, kept, context)?;
    for (column, value) in values.iter().enumerate() {
        let scale = (divisor * value).sqrt().recip();
        for row in 0..matrix.cols() {
            axes[(row, column)] *= scale;
        }
    }
    Ok((values, axes))
}

fn canonicalize_signs(components: &mut Mat<f64>) {
    for column in 0..components.ncols() {
        let mut pivot: Option<(usize, f64)> = None;
        for row in 0..components.nrows() {
            let value = components[(row, column)];
            // Largest magnitude wins; the earliest index wins ties.
            if pivot.is_none_or(|(_, best)| value.abs().total_cmp(&best.abs()).is_gt()) {
                pivot = Some((row, value));
            }
        }
        if pivot.is_some_and(|(_, value)| value < 0.0) {
            for row in 0..components.nrows() {
                components[(row, column)] = -components[(row, column)];
            }
        }
    }
}

#[cfg(test)]
#[path = "pca_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "pca_reference_tests.rs"]
mod reference_tests;
