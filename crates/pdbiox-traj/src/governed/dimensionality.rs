//! Governed dimensionality reduction.

use pdbiox_core::contract::{AlgorithmId, Analysis, AnalysisPolicy, Coverage, ParameterValue};
use pdbiox_geom::PeriodicAngle;

use crate::{
    CartesianFit, DiffusionMap, EnsembleDistanceMatrix, EnsembleGeometryError, FrameView,
    PcaResult, cartesian_pca, cartesian_pca_view, diffusion_map, dihedral_pca,
};

/// Runs Cartesian PCA and records every numerical choice in provenance.
///
/// # Errors
///
/// Returns the underlying validation, fitting, decomposition or memory error.
pub fn analyse_cartesian_pca(
    frames: &[Vec<[f32; 3]>],
    fit: CartesianFit<'_>,
    components: usize,
    memory_limit: usize,
    policy: &AnalysisPolicy,
) -> Result<Analysis<PcaResult>, EnsembleGeometryError> {
    cartesian_pca_analysis(
        cartesian_pca(frames, fit, components, memory_limit)?,
        frames.len(),
        fit,
        components,
        memory_limit,
        policy,
    )
}

/// Runs Cartesian PCA over borrowed contiguous frames with governed provenance.
///
/// # Errors
///
/// Returns the underlying validation, fitting, decomposition or memory error.
pub fn analyse_cartesian_pca_view(
    frames: FrameView<'_>,
    fit: CartesianFit<'_>,
    components: usize,
    memory_limit: usize,
    policy: &AnalysisPolicy,
) -> Result<Analysis<PcaResult>, EnsembleGeometryError> {
    cartesian_pca_analysis(
        cartesian_pca_view(frames, fit, components, memory_limit)?,
        frames.frame_count(),
        fit,
        components,
        memory_limit,
        policy,
    )
}

fn cartesian_pca_analysis(
    value: PcaResult,
    frame_count: usize,
    fit: CartesianFit<'_>,
    components: usize,
    memory_limit: usize,
    policy: &AnalysisPolicy,
) -> Result<Analysis<PcaResult>, EnsembleGeometryError> {
    let fit_name = match fit {
        CartesianFit::None => "none",
        CartesianFit::Reference(_) => "reference",
        CartesianFit::IterativeMean { .. } => "iterative-mean",
    };
    let mut analysis = Analysis::complete(value, coverage(frame_count)?, policy);
    analysis.provenance = analysis
        .provenance
        .with_algorithm(AlgorithmId::new("cartesian-pca", "1"))
        .with_parameter("fit", ParameterValue::Text(fit_name.into()))
        .with_parameter("components", integer(components)?)
        .with_parameter("memory_limit", integer(memory_limit)?);
    Ok(analysis)
}

/// Runs circular-feature dihedral PCA with governed provenance.
///
/// # Errors
///
/// Returns the underlying dimensionality, decomposition or memory error.
pub fn analyse_dihedral_pca(
    observations: &[Vec<PeriodicAngle>],
    torsion_set: &str,
    components: usize,
    memory_limit: usize,
    policy: &AnalysisPolicy,
) -> Result<Analysis<PcaResult>, EnsembleGeometryError> {
    let value = dihedral_pca(observations, components, memory_limit)?;
    let mut analysis = Analysis::complete(value, coverage(observations.len())?, policy);
    analysis.provenance = analysis
        .provenance
        .with_algorithm(AlgorithmId::new("dihedral-pca-cos-sin", "1"))
        .with_parameter("torsion_set", ParameterValue::Text(torsion_set.into()))
        .with_parameter("components", integer(components)?)
        .with_parameter("memory_limit", integer(memory_limit)?);
    Ok(analysis)
}

/// Runs a governed diffusion map over an explicitly named metric matrix.
///
/// # Errors
///
/// Returns the underlying matrix or diffusion-kernel validation error.
pub fn analyse_diffusion_map(
    distances: &EnsembleDistanceMatrix,
    metric: &str,
    epsilon: f64,
    time: u32,
    dimensions: usize,
    policy: &AnalysisPolicy,
) -> Result<Analysis<DiffusionMap>, EnsembleGeometryError> {
    let value = diffusion_map(distances, epsilon, time, dimensions)?;
    let mut analysis = Analysis::complete(value, coverage(distances.size)?, policy);
    let Some(epsilon) = ParameterValue::finite_float(epsilon) else {
        return Err(EnsembleGeometryError::InvalidParameter);
    };
    analysis.provenance = analysis
        .provenance
        .with_algorithm(AlgorithmId::new("diffusion-map-row-normalized", "1"))
        .with_parameter("metric", ParameterValue::Text(metric.into()))
        .with_parameter("epsilon", epsilon)
        .with_parameter("time", ParameterValue::Integer(i64::from(time)))
        .with_parameter("dimensions", integer(dimensions)?);
    Ok(analysis)
}

fn coverage(observations: usize) -> Result<Coverage, EnsembleGeometryError> {
    u32::try_from(observations)
        .map(Coverage::complete)
        .map_err(|_| EnsembleGeometryError::InvalidParameter)
}

fn integer(value: usize) -> Result<ParameterValue, EnsembleGeometryError> {
    i64::try_from(value)
        .map(ParameterValue::Integer)
        .map_err(|_| EnsembleGeometryError::InvalidParameter)
}
