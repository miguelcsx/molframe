//! Direct Python entry points for trajectory kernels that are not readers.

use crate::contract::{PyAnalysis, analysis_with_value};
use crate::intrinsic::{
    PyCartesianFit, PyClustering, PyConvergenceBlock, PyEnsembleDistanceMatrix, PyFrameAlignment,
    PyGroupVariance, PyHarmonicSimilarity, PyHarmonicSimilarityOptions, PyLinkage,
    PyRemainderPolicy, borrowed_frames,
};
use crate::query::PyAnalysisPolicy;
use numpy::{PyReadonlyArray2, PyReadonlyArray3};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyList;
use std::fmt;

#[derive(Debug)]
enum NativeFrameAnalysisError {
    Trajectory(pdbiox::traj::TrajectoryError),
    Geometry(String),
    Invalid(String),
}

impl From<pdbiox::traj::TrajectoryError> for NativeFrameAnalysisError {
    fn from(value: pdbiox::traj::TrajectoryError) -> Self {
        Self::Trajectory(value)
    }
}

impl fmt::Display for NativeFrameAnalysisError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Trajectory(error) => error.fmt(formatter),
            Self::Geometry(error) | Self::Invalid(error) => formatter.write_str(error),
        }
    }
}

struct NativeRmsdAnalysis {
    reference: usize,
    alignment: pdbiox::traj::FrameAlignment,
    policy: pdbiox::AnalysisPolicy,
}

struct RmsdPartial {
    reference: Box<[[f32; 3]]>,
    values: Vec<f64>,
}

impl pdbiox::traj::FrameAnalysis for NativeRmsdAnalysis {
    type Output = Vec<f64>;
    type Partial = RmsdPartial;
    type Error = NativeFrameAnalysisError;

    const NAME: &'static str = "rmsd_to_reference";
    const PARALLELIZABLE: bool = true;

    fn prepare(&self, trajectory: &pdbiox::traj::Trajectory) -> Result<Self::Partial, Self::Error> {
        let frame = trajectory
            .frame(self.reference)
            .ok_or_else(|| NativeFrameAnalysisError::Invalid("reference frame is absent".into()))?;
        Ok(RmsdPartial {
            reference: frame.positions.to_vec().into_boxed_slice(),
            values: Vec::new(),
        })
    }

    fn single_frame(
        &self,
        timestep: &pdbiox::traj::Timestep,
        partial: &mut Self::Partial,
        _context: &pdbiox::core::ExecutionContext,
    ) -> Result<(), Self::Error> {
        if timestep.positions.len() != partial.reference.len() {
            return Err(pdbiox::traj::TrajectoryError::AtomCountMismatch {
                expected: partial.reference.len(),
                found: timestep.positions.len(),
            }
            .into());
        }
        let value = match self.alignment {
            pdbiox::traj::FrameAlignment::None => {
                { pdbiox::rmsd(&timestep.positions, &partial.reference) }
                    .map_err(|error| NativeFrameAnalysisError::Geometry(format!("{error:?}")))?
            }
            pdbiox::traj::FrameAlignment::Rigid => {
                pdbiox::superpose(&timestep.positions, &partial.reference)
                    .map(|fit| fit.rmsd)
                    .map_err(|error| NativeFrameAnalysisError::Geometry(format!("{error:?}")))?
            }
        };
        partial.values.push(value);
        Ok(())
    }

    fn merge(&self, partials: Vec<Self::Partial>) -> Result<Self::Partial, Self::Error> {
        let mut partials = partials.into_iter();
        let Some(mut merged) = partials.next() else {
            return Err(NativeFrameAnalysisError::Invalid(
                "frame analysis produced no blocks".into(),
            ));
        };
        for partial in partials {
            merged.values.extend(partial.values);
        }
        Ok(merged)
    }

    fn conclude(
        &self,
        partial: Self::Partial,
    ) -> Result<pdbiox::Analysis<Self::Output>, Self::Error> {
        let count = u32::try_from(partial.values.len()).map_err(|_| {
            NativeFrameAnalysisError::Invalid("frame count exceeds analysis coverage limits".into())
        })?;
        let mut result = pdbiox::Analysis::complete(
            partial.values,
            pdbiox::Coverage::complete(count),
            &self.policy,
        );
        result.provenance = result
            .provenance
            .with_algorithm(pdbiox::AlgorithmId::new(Self::NAME, "1"))
            .with_parameter(
                "reference",
                pdbiox::ParameterValue::Integer(i64::try_from(self.reference).map_err(|_| {
                    NativeFrameAnalysisError::Invalid("reference exceeds parameter limits".into())
                })?),
            )
            .with_parameter(
                "alignment",
                pdbiox::ParameterValue::Text(format!("{:?}", self.alignment).into()),
            );
        Ok(result)
    }
}

#[derive(Clone, Copy, Debug)]
enum FrameAnalysisKind {
    RmsdToReference {
        reference: usize,
        alignment: pdbiox::traj::FrameAlignment,
    },
}

#[pyclass(name = "FrameAnalysis", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyFrameAnalysis {
    kind: FrameAnalysisKind,
}

#[pymethods]
impl PyFrameAnalysis {
    #[staticmethod]
    fn rmsd_to_reference(reference: usize, alignment: PyFrameAlignment) -> Self {
        Self {
            kind: FrameAnalysisKind::RmsdToReference {
                reference,
                alignment: alignment.into(),
            },
        }
    }

    #[getter]
    fn name(&self) -> &'static str {
        match self.kind {
            FrameAnalysisKind::RmsdToReference { .. } => "rmsd_to_reference",
        }
    }
}

impl PyFrameAnalysis {
    fn native(&self, policy: pdbiox::AnalysisPolicy) -> NativeRmsdAnalysis {
        match self.kind {
            FrameAnalysisKind::RmsdToReference {
                reference,
                alignment,
            } => NativeRmsdAnalysis {
                reference,
                alignment,
                policy,
            },
        }
    }
}

#[pyfunction]
#[pyo3(signature = (trajectory, analysis, *, policy, context=None))]
pub(crate) fn run_analysis(
    py: Python<'_>,
    trajectory: PyRef<'_, super::PyTrajectory>,
    analysis: &PyFrameAnalysis,
    policy: &PyAnalysisPolicy,
    context: Option<&crate::core::execution::PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let native_trajectory = trajectory.native_trajectory(py)?;
    let native_analysis = analysis.native(policy.inner.clone());
    let context = context.map_or_else(
        crate::core::execution::default_context,
        crate::core::execution::PyExecutionContext::native,
    );
    let result = py
        .detach(move || pdbiox::traj::run_analysis(&native_trajectory, &native_analysis, &context))
        .map_err(kernel_error)?;
    analysis_with_value(py, result, |py, value| {
        Ok(PyList::new(py, value)?.unbind().into_any())
    })
}

#[pyclass(name = "PathFrameMetric", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyPathFrameMetric {
    CartesianRmsd,
    FittedRmsd,
}

#[pyclass(name = "PathSimilarity", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyPathSimilarity {
    #[pyo3(get)]
    pub(crate) hausdorff_distance: f64,
    #[pyo3(get)]
    pub(crate) discrete_frechet_distance: f64,
}

#[pyclass(name = "SurvivalMode", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PySurvivalMode {
    Continuous,
    Intermittent,
}

#[pyclass(name = "WaterSurvival", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyWaterSurvival {
    #[pyo3(get)]
    pub(crate) lag: usize,
    #[pyo3(get)]
    pub(crate) origins: u64,
    #[pyo3(get)]
    pub(crate) probability: Option<f64>,
}

#[pyclass(name = "WaterDynamics", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyWaterDynamics {
    #[pyo3(get)]
    pub(crate) survival: Vec<PyWaterSurvival>,
    #[pyo3(get)]
    pub(crate) residence_time: Option<f64>,
}

#[pyfunction]
pub(crate) fn agglomerative_clustering(
    py: Python<'_>,
    distances: &PyEnsembleDistanceMatrix,
    cluster_count: usize,
    linkage: PyLinkage,
) -> PyResult<PyClustering> {
    py.detach(move || -> PyResult<PyClustering> {
        pdbiox::traj::agglomerative_clustering(&distances.native(), cluster_count, linkage.into())
            .map(PyClustering::from)
            .map_err(kernel_error)
    })
}

#[pyfunction]
pub(crate) fn dbscan_clustering(
    py: Python<'_>,
    distances: &PyEnsembleDistanceMatrix,
    epsilon: f64,
    minimum_points: usize,
) -> PyResult<PyClustering> {
    py.detach(move || -> PyResult<PyClustering> {
        pdbiox::traj::dbscan_clustering(&distances.native(), epsilon, minimum_points)
            .map(PyClustering::from)
            .map_err(kernel_error)
    })
}

#[pyfunction]
pub(crate) fn medoid(
    py: Python<'_>,
    distances: &PyEnsembleDistanceMatrix,
    members: Vec<usize>,
) -> PyResult<usize> {
    py.detach(move || -> PyResult<usize> {
        pdbiox::traj::medoid(&distances.native(), &members).map_err(kernel_error)
    })
}

#[pyfunction]
pub(crate) fn harmonic_ensemble_similarity(
    py: Python<'_>,
    first: Vec<Vec<f64>>,
    second: Vec<Vec<f64>>,
    options: &PyHarmonicSimilarityOptions,
) -> PyResult<PyHarmonicSimilarity> {
    py.detach(move || -> PyResult<PyHarmonicSimilarity> {
        pdbiox::traj::harmonic_ensemble_similarity(&first, &second, options.native())
            .map(PyHarmonicSimilarity::from)
            .map_err(kernel_error)
    })
}

#[pyfunction]
pub(crate) fn cluster_population_similarity(
    py: Python<'_>,
    first: Vec<usize>,
    second: Vec<usize>,
    cluster_count: usize,
) -> PyResult<f64> {
    py.detach(move || -> PyResult<f64> {
        pdbiox::traj::cluster_population_similarity(&first, &second, cluster_count)
            .map_err(kernel_error)
    })
}

#[pyfunction]
pub(crate) fn group_coordinate_variance(
    py: Python<'_>,
    frames: PyReadonlyArray3<'_, f32>,
    groups: Vec<Vec<usize>>,
) -> PyResult<Vec<PyGroupVariance>> {
    let frames = borrowed_frames(&frames)?;
    py.detach(move || pdbiox::traj::group_coordinate_variance_view(frames, &groups))
        .map(|values| values.into_iter().map(PyGroupVariance::from).collect())
        .map_err(kernel_error)
}

#[pyfunction]
pub(crate) fn block_convergence(
    py: Python<'_>,
    values: Vec<f64>,
    block_size: usize,
    remainder: PyRemainderPolicy,
) -> PyResult<Vec<PyConvergenceBlock>> {
    py.detach(move || -> PyResult<Vec<PyConvergenceBlock>> {
        pdbiox::traj::block_convergence(&values, block_size, remainder.into())
            .map(|values| values.into_iter().map(PyConvergenceBlock::from).collect())
            .map_err(kernel_error)
    })
}

#[pyfunction]
pub(crate) fn path_similarity(
    py: Python<'_>,
    first: PyReadonlyArray3<'_, f32>,
    second: PyReadonlyArray3<'_, f32>,
    metric: PyPathFrameMetric,
    memory_limit_bytes: usize,
) -> PyResult<PyPathSimilarity> {
    let first = borrowed_frames(&first)?;
    let second = borrowed_frames(&second)?;
    let metric = match metric {
        PyPathFrameMetric::CartesianRmsd => pdbiox::traj::PathFrameMetric::CartesianRmsd,
        PyPathFrameMetric::FittedRmsd => pdbiox::traj::PathFrameMetric::FittedRmsd,
    };
    py.detach(move || pdbiox::traj::path_similarity_view(first, second, metric, memory_limit_bytes))
        .map(|value| PyPathSimilarity {
            hausdorff_distance: value.hausdorff_distance,
            discrete_frechet_distance: value.discrete_frechet_distance,
        })
        .map_err(kernel_error)
}

#[pyfunction(name = "trajectory_water_dynamics")]
pub(crate) fn water_dynamics(
    occupancy: Vec<Vec<bool>>,
    maximum_lag: usize,
    frame_duration: f64,
    mode: PySurvivalMode,
) -> PyResult<PyWaterDynamics> {
    let mode = match mode {
        PySurvivalMode::Continuous => pdbiox::traj::SurvivalMode::Continuous,
        PySurvivalMode::Intermittent => pdbiox::traj::SurvivalMode::Intermittent,
    };
    pdbiox::traj::water_dynamics(&occupancy, maximum_lag, frame_duration, mode)
        .map(|value| PyWaterDynamics {
            survival: value
                .survival
                .into_iter()
                .map(|entry| PyWaterSurvival {
                    lag: entry.lag,
                    origins: entry.origins,
                    probability: entry.probability,
                })
                .collect(),
            residence_time: value.residence_time,
        })
        .map_err(kernel_error)
}

#[pyfunction]
#[pyo3(signature = (frames, *, fit, components, memory_limit, policy))]
pub(crate) fn analyse_cartesian_pca(
    py: Python<'_>,
    frames: PyReadonlyArray3<'_, f32>,
    fit: &PyCartesianFit,
    components: usize,
    memory_limit: usize,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    crate::intrinsic::analyse_pca(py, frames, fit, components, memory_limit, policy)
}

#[pyfunction]
#[pyo3(signature = (angles, *, torsion_set, components, memory_limit, policy))]
pub(crate) fn analyse_dihedral_pca(
    py: Python<'_>,
    angles: PyReadonlyArray2<'_, f64>,
    torsion_set: &str,
    components: usize,
    memory_limit: usize,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    crate::intrinsic::analyse_torsion_pca(py, angles, torsion_set, components, memory_limit, policy)
}

#[pyfunction]
#[pyo3(signature = (distances, *, metric, epsilon, time, dimensions, policy))]
pub(crate) fn analyse_diffusion_map(
    py: Python<'_>,
    distances: PyReadonlyArray2<'_, f64>,
    metric: &str,
    epsilon: f64,
    time: u32,
    dimensions: usize,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    crate::intrinsic::analyse_diffusion(py, distances, metric, epsilon, time, dimensions, policy)
}

fn kernel_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add(
        "DEFAULT_PAIRWISE_MEMORY_LIMIT",
        pdbiox::traj::DEFAULT_PAIRWISE_MEMORY_LIMIT,
    )?;
    module.add(
        "DEFAULT_PROCRUSTES_TOLERANCE",
        pdbiox::traj::DEFAULT_PROCRUSTES_TOLERANCE,
    )?;
    module.add(
        "DEFAULT_PROCRUSTES_MAXIMUM_ITERATIONS",
        pdbiox::traj::DEFAULT_PROCRUSTES_MAXIMUM_ITERATIONS,
    )?;
    module.add_class::<PyFrameAnalysis>()?;
    module.add_class::<PyPathFrameMetric>()?;
    module.add_class::<PyPathSimilarity>()?;
    module.add_class::<PySurvivalMode>()?;
    module.add_class::<PyWaterSurvival>()?;
    module.add_class::<PyWaterDynamics>()?;
    module.add_function(wrap_pyfunction!(agglomerative_clustering, module)?)?;
    module.add_function(wrap_pyfunction!(dbscan_clustering, module)?)?;
    module.add_function(wrap_pyfunction!(medoid, module)?)?;
    module.add_function(wrap_pyfunction!(harmonic_ensemble_similarity, module)?)?;
    module.add_function(wrap_pyfunction!(cluster_population_similarity, module)?)?;
    module.add_function(wrap_pyfunction!(group_coordinate_variance, module)?)?;
    module.add_function(wrap_pyfunction!(block_convergence, module)?)?;
    module.add_function(wrap_pyfunction!(path_similarity, module)?)?;
    module.add_function(wrap_pyfunction!(water_dynamics, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_cartesian_pca, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_dihedral_pca, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_diffusion_map, module)?)?;
    module.add_function(wrap_pyfunction!(run_analysis, module)?)?;
    Ok(())
}
