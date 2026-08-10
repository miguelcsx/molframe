//! Governed wrappers for deterministic ensemble kernels.

use crate::{
    Clustering, ConvergenceBlock, EnsembleDistanceMatrix, EnsembleGeometryError,
    EnsembleSimilarityError, EnsembleStatisticsError, FrameAlignment, GroupVariance,
    HarmonicSimilarity, HarmonicSimilarityOptions, KMeans, KMeansError, KMeansOptions, Linkage,
    RemainderPolicy, Timestep, agglomerative_clustering, block_convergence,
    cluster_population_similarity, dbscan_clustering, generalized_procrustes_mean,
    group_coordinate_variance, harmonic_ensemble_similarity, kmeans, medoid, pairwise_fitted_rmsd,
    pairwise_torus_distance, rmsd_to_reference,
};
use pdbiox_core::contract::{AlgorithmId, Analysis, AnalysisPolicy, Coverage, ParameterValue};
use pdbiox_geom::{PeriodicAngle, TorusMetric};

/// Raw ensemble-domain or governance failure.
#[derive(Debug, thiserror::Error)]
pub enum GovernedEnsembleError {
    /// Pairwise geometry, clustering or representative failure.
    #[error(transparent)]
    Geometry(#[from] EnsembleGeometryError),
    /// Harmonic or population-similarity failure.
    #[error(transparent)]
    Similarity(#[from] EnsembleSimilarityError),
    /// K-means failure.
    #[error(transparent)]
    KMeans(#[from] KMeansError),
    /// Variance or convergence failure.
    #[error(transparent)]
    Statistics(#[from] EnsembleStatisticsError),
    /// Observation count exceeds public coverage counters.
    #[error("ensemble observation count exceeds u32 coverage limits")]
    CoverageOverflow,
}

fn governed<T>(
    value: T,
    observations: usize,
    policy: &AnalysisPolicy,
    algorithm: &'static str,
) -> Result<Analysis<T>, GovernedEnsembleError> {
    let count = u32::try_from(observations).map_err(|_| GovernedEnsembleError::CoverageOverflow)?;
    let mut analysis = Analysis::complete(value, Coverage::complete(count), policy);
    analysis.provenance = analysis
        .provenance
        .with_algorithm(AlgorithmId::new(algorithm, "1"));
    Ok(analysis)
}

fn integer(value: usize) -> Result<ParameterValue, GovernedEnsembleError> {
    i64::try_from(value)
        .map(ParameterValue::Integer)
        .map_err(|_| GovernedEnsembleError::CoverageOverflow)
}

fn float(value: f64) -> Result<ParameterValue, GovernedEnsembleError> {
    ParameterValue::finite_float(value).ok_or(GovernedEnsembleError::Geometry(
        EnsembleGeometryError::InvalidParameter,
    ))
}

/// Governed RMSD series against an explicit reference frame.
///
/// # Errors
///
/// Returns the typed raw-kernel failure or coverage overflow.
pub fn analyse_rmsd_to_reference(
    frames: &[Timestep],
    reference: usize,
    alignment: FrameAlignment,
    policy: &AnalysisPolicy,
) -> Result<Analysis<Vec<f64>>, GovernedEnsembleError> {
    let value = rmsd_to_reference(frames, reference, alignment)?;
    let mut result = governed(value, frames.len(), policy, "ensemble-rmsd-to-reference")?;
    result.provenance = result
        .provenance
        .with_parameter("reference", integer(reference)?)
        .with_parameter(
            "alignment",
            ParameterValue::Text(format!("{alignment:?}").into()),
        );
    Ok(result)
}

/// Governed pairwise fitted RMSD matrix.
///
/// # Errors
///
/// Returns the typed raw-kernel failure or coverage overflow.
pub fn analyse_pairwise_fitted_rmsd(
    frames: &[Vec<[f32; 3]>],
    memory_limit: usize,
    policy: &AnalysisPolicy,
) -> Result<Analysis<EnsembleDistanceMatrix>, GovernedEnsembleError> {
    let value = pairwise_fitted_rmsd(frames, memory_limit)?;
    let mut result = governed(value, frames.len(), policy, "pairwise-fitted-rmsd")?;
    result.provenance = result
        .provenance
        .with_parameter("memory_limit", integer(memory_limit)?);
    Ok(result)
}

/// Governed pairwise torus-distance matrix.
///
/// # Errors
///
/// Returns the typed raw-kernel failure or coverage overflow.
pub fn analyse_pairwise_torus_distance(
    points: &[Vec<PeriodicAngle>],
    metric: &TorusMetric,
    metric_name: &str,
    memory_limit: usize,
    policy: &AnalysisPolicy,
) -> Result<Analysis<EnsembleDistanceMatrix>, GovernedEnsembleError> {
    let value = pairwise_torus_distance(points, metric, memory_limit)?;
    let mut result = governed(value, points.len(), policy, "pairwise-torus-distance")?;
    result.provenance = result
        .provenance
        .with_parameter("metric", ParameterValue::Text(metric_name.into()))
        .with_parameter("memory_limit", integer(memory_limit)?);
    Ok(result)
}

/// Governed generalized Procrustes centroid.
///
/// # Errors
///
/// Returns the typed raw-kernel failure or coverage overflow.
pub fn analyse_generalized_procrustes_mean(
    frames: &[Vec<[f32; 3]>],
    tolerance: f64,
    maximum_iterations: usize,
    policy: &AnalysisPolicy,
) -> Result<Analysis<Vec<[f32; 3]>>, GovernedEnsembleError> {
    let value = generalized_procrustes_mean(frames, tolerance, maximum_iterations)?;
    let mut result = governed(value, frames.len(), policy, "generalized-procrustes-mean")?;
    result.provenance = result
        .provenance
        .with_parameter("tolerance", float(tolerance)?)
        .with_parameter("maximum_iterations", integer(maximum_iterations)?);
    Ok(result)
}

/// Governed agglomerative clustering.
///
/// # Errors
///
/// Returns the typed raw-kernel failure or coverage overflow.
pub fn analyse_agglomerative_clustering(
    distances: &EnsembleDistanceMatrix,
    cluster_count: usize,
    linkage: Linkage,
    policy: &AnalysisPolicy,
) -> Result<Analysis<Clustering>, GovernedEnsembleError> {
    let value = agglomerative_clustering(distances, cluster_count, linkage)?;
    let mut result = governed(value, distances.size, policy, "agglomerative-clustering")?;
    result.provenance = result
        .provenance
        .with_parameter("cluster_count", integer(cluster_count)?)
        .with_parameter(
            "linkage",
            ParameterValue::Text(format!("{linkage:?}").into()),
        );
    Ok(result)
}

/// Governed DBSCAN clustering.
///
/// # Errors
///
/// Returns the typed raw-kernel failure or coverage overflow.
pub fn analyse_dbscan_clustering(
    distances: &EnsembleDistanceMatrix,
    epsilon: f64,
    minimum_points: usize,
    policy: &AnalysisPolicy,
) -> Result<Analysis<Clustering>, GovernedEnsembleError> {
    let value = dbscan_clustering(distances, epsilon, minimum_points)?;
    let mut result = governed(value, distances.size, policy, "dbscan-clustering")?;
    result.provenance = result
        .provenance
        .with_parameter("epsilon", float(epsilon)?)
        .with_parameter("minimum_points", integer(minimum_points)?);
    Ok(result)
}

/// Governed medoid representative selection.
///
/// # Errors
///
/// Returns the typed raw-kernel failure or coverage overflow.
pub fn analyse_medoid(
    distances: &EnsembleDistanceMatrix,
    members: &[usize],
    policy: &AnalysisPolicy,
) -> Result<Analysis<usize>, GovernedEnsembleError> {
    governed(
        medoid(distances, members)?,
        members.len(),
        policy,
        "ensemble-medoid",
    )
}

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
    let value = kmeans(observations, options)?;
    let mut result = governed(value, observations.len(), policy, "kmeans")?;
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

/// Governed harmonic ensemble similarity.
///
/// # Errors
///
/// Returns the typed raw-kernel failure or coverage overflow.
pub fn analyse_harmonic_ensemble_similarity(
    first: &[Vec<f64>],
    second: &[Vec<f64>],
    options: HarmonicSimilarityOptions,
    policy: &AnalysisPolicy,
) -> Result<Analysis<HarmonicSimilarity>, GovernedEnsembleError> {
    let value = harmonic_ensemble_similarity(first, second, options)?;
    let mut result = governed(
        value,
        first.len() + second.len(),
        policy,
        "harmonic-ensemble-similarity",
    )?;
    result.provenance = result
        .provenance
        .with_parameter(
            "covariance_regularization",
            float(options.covariance_regularization)?,
        )
        .with_parameter("memory_limit_bytes", integer(options.memory_limit_bytes)?);
    Ok(result)
}

/// Governed cluster-population similarity.
///
/// # Errors
///
/// Returns the typed raw-kernel failure or coverage overflow.
pub fn analyse_cluster_population_similarity(
    first: &[usize],
    second: &[usize],
    cluster_count: usize,
    policy: &AnalysisPolicy,
) -> Result<Analysis<f64>, GovernedEnsembleError> {
    let value = cluster_population_similarity(first, second, cluster_count)?;
    let mut result = governed(
        value,
        first.len() + second.len(),
        policy,
        "cluster-population-similarity",
    )?;
    result.provenance = result
        .provenance
        .with_parameter("cluster_count", integer(cluster_count)?);
    Ok(result)
}

/// Governed per-group coordinate variance.
///
/// # Errors
///
/// Returns the typed raw-kernel failure or coverage overflow.
pub fn analyse_group_coordinate_variance(
    frames: &[Vec<[f32; 3]>],
    groups: &[Vec<usize>],
    policy: &AnalysisPolicy,
) -> Result<Analysis<Vec<GroupVariance>>, GovernedEnsembleError> {
    governed(
        group_coordinate_variance(frames, groups)?,
        frames.len(),
        policy,
        "group-coordinate-variance",
    )
}

/// Governed scalar block-convergence summary.
///
/// # Errors
///
/// Returns the typed raw-kernel failure or coverage overflow.
pub fn analyse_block_convergence(
    values: &[f64],
    block_size: usize,
    remainder: RemainderPolicy,
    policy: &AnalysisPolicy,
) -> Result<Analysis<Vec<ConvergenceBlock>>, GovernedEnsembleError> {
    let value = block_convergence(values, block_size, remainder)?;
    let mut result = governed(value, values.len(), policy, "block-convergence")?;
    result.provenance = result
        .provenance
        .with_parameter("block_size", integer(block_size)?)
        .with_parameter(
            "remainder",
            ParameterValue::Text(format!("{remainder:?}").into()),
        );
    Ok(result)
}

#[cfg(test)]
#[path = "ensemble_tests.rs"]
mod tests;
