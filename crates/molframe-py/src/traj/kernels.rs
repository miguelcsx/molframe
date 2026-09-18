//! Direct trajectory kernels over contiguous `NumPy` arrays.

use crate::intrinsic::{
    PyEnsembleDistanceMatrix, PyFrameAlignment, PyKMeans, PyKMeansOptions, borrowed_frames,
};
use numpy::ndarray::Array2;
use numpy::{
    IntoPyArray, PyArray2, PyReadonlyArray1, PyReadonlyArray2, PyReadonlyArray3,
    PyUntypedArrayMethods,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

const DEFAULT_MEMORY_LIMIT: usize = 512 * 1024 * 1024;
const CONTIGUOUS: &str = "trajectory arrays must be C-contiguous; pass an explicit contiguous copy";

#[pyclass(name = "MeanSquaredDisplacement", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyMeanSquaredDisplacement {
    #[pyo3(get)]
    pub(crate) lag: usize,
    #[pyo3(get)]
    pub(crate) observations: u64,
    #[pyo3(get)]
    pub(crate) value: f64,
}

impl From<molframe::traj::MeanSquaredDisplacement> for PyMeanSquaredDisplacement {
    fn from(value: molframe::traj::MeanSquaredDisplacement) -> Self {
        Self {
            lag: value.lag,
            observations: value.observations,
            value: value.value,
        }
    }
}

#[pyclass(name = "TrajectoryDielectricOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyTrajectoryDielectricOptions {
    #[pyo3(get)]
    pub(crate) fluctuation_prefactor: f64,
}

#[pymethods]
impl PyTrajectoryDielectricOptions {
    #[new]
    fn new(fluctuation_prefactor: f64) -> PyResult<Self> {
        if !fluctuation_prefactor.is_finite() || fluctuation_prefactor < 0.0 {
            return Err(PyValueError::new_err(
                "fluctuation_prefactor must be finite and non-negative",
            ));
        }
        Ok(Self {
            fluctuation_prefactor,
        })
    }
}

#[pyclass(name = "DielectricEstimate", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyDielectricEstimate {
    #[pyo3(get)]
    pub(crate) mean_dipole: [f64; 3],
    #[pyo3(get)]
    pub(crate) covariance: [[f64; 3]; 3],
    #[pyo3(get)]
    pub(crate) fluctuation: f64,
    #[pyo3(get)]
    pub(crate) relative_permittivity: f64,
}

impl From<molframe::traj::DielectricEstimate> for PyDielectricEstimate {
    fn from(value: molframe::traj::DielectricEstimate) -> Self {
        Self {
            mean_dipole: value.mean_dipole,
            covariance: value.covariance,
            fluctuation: value.fluctuation,
            relative_permittivity: value.relative_permittivity,
        }
    }
}

#[pyfunction]
#[pyo3(signature = (observations, options))]
pub(crate) fn kmeans(
    py: Python<'_>,
    observations: PyReadonlyArray2<'_, f64>,
    options: &PyKMeansOptions,
) -> PyResult<PyKMeans> {
    let shape = observations.shape();
    let rows = shape[0];
    let columns = shape[1];
    let observations = observations
        .as_slice()
        .map_err(|_| PyValueError::new_err(CONTIGUOUS))?;
    let options = options.clone();
    py.detach(move || {
        molframe::traj::kmeans_view(observations, rows, columns, options.native())
            .map(PyKMeans::from)
            .map_err(value_error)
    })
}

#[pyfunction]
#[pyo3(signature = (frames, maximum_lag, *, atoms=None))]
pub(crate) fn mean_squared_displacement(
    py: Python<'_>,
    frames: PyReadonlyArray3<'_, f32>,
    maximum_lag: usize,
    atoms: Option<PyReadonlyArray1<'_, usize>>,
) -> PyResult<Vec<PyMeanSquaredDisplacement>> {
    let view = borrowed_frames(&frames)?;
    let atoms = match atoms.as_ref() {
        Some(values) => values
            .as_slice()
            .map_err(|_| PyValueError::new_err(CONTIGUOUS))?,
        None => &[],
    };
    py.detach(move || {
        molframe::traj::mean_squared_displacement_view(view, atoms, maximum_lag)
            .map(|values| values.into_iter().map(Into::into).collect())
            .map_err(value_error)
    })
}

#[pyfunction(name = "trajectory_dielectric_from_dipoles")]
#[pyo3(signature = (dipoles, options))]
pub(crate) fn trajectory_dielectric_from_dipoles(
    py: Python<'_>,
    dipoles: PyReadonlyArray2<'_, f64>,
    options: PyTrajectoryDielectricOptions,
) -> PyResult<PyDielectricEstimate> {
    let dipoles = triplets(dipoles, "dipoles must have shape (observations, 3)")?;
    py.detach(move || {
        molframe::traj::dielectric_from_dipoles(
            &dipoles,
            molframe::traj::DielectricOptions {
                fluctuation_prefactor: options.fluctuation_prefactor,
            },
        )
        .map(Into::into)
        .map_err(value_error)
    })
}

#[pyfunction]
#[pyo3(signature = (frames, reference, alignment))]
pub(crate) fn rmsd_to_reference(
    py: Python<'_>,
    frames: PyReadonlyArray3<'_, f32>,
    reference: usize,
    alignment: PyFrameAlignment,
) -> PyResult<Vec<f64>> {
    let view = borrowed_frames(&frames)?;
    py.detach(move || {
        molframe::traj::rmsd_to_reference_view(view, reference, alignment.into())
            .map_err(value_error)
    })
}

#[pyfunction]
#[pyo3(signature = (frames, memory_limit=DEFAULT_MEMORY_LIMIT))]
pub(crate) fn pairwise_fitted_rmsd(
    py: Python<'_>,
    frames: PyReadonlyArray3<'_, f32>,
    memory_limit: usize,
) -> PyResult<PyEnsembleDistanceMatrix> {
    let view = borrowed_frames(&frames)?;
    py.detach(move || {
        molframe::traj::pairwise_fitted_rmsd_view(view, memory_limit)
            .map(Into::into)
            .map_err(value_error)
    })
}

#[pyfunction]
#[pyo3(signature = (points, weights, memory_limit=DEFAULT_MEMORY_LIMIT))]
pub(crate) fn pairwise_torus_distance(
    py: Python<'_>,
    points: PyReadonlyArray2<'_, f64>,
    weights: PyReadonlyArray1<'_, f64>,
    memory_limit: usize,
) -> PyResult<PyEnsembleDistanceMatrix> {
    let points = periodic_points(points)?;
    let weights = weights
        .as_slice()
        .map_err(|_| PyValueError::new_err(CONTIGUOUS))?
        .to_vec();
    let metric = molframe::TorusMetric::new(weights).map_err(value_error)?;
    py.detach(move || {
        molframe::traj::pairwise_torus_distance(&points, &metric, memory_limit)
            .map(Into::into)
            .map_err(value_error)
    })
}

#[pyfunction]
#[pyo3(signature = (frames, tolerance, maximum_iterations))]
pub(crate) fn generalized_procrustes_mean(
    py: Python<'_>,
    frames: PyReadonlyArray3<'_, f32>,
    tolerance: f64,
    maximum_iterations: usize,
) -> PyResult<Py<PyArray2<f32>>> {
    let view = borrowed_frames(&frames)?;
    let result = py
        .detach(move || {
            molframe::traj::generalized_procrustes_mean_view(view, tolerance, maximum_iterations)
        })
        .map_err(value_error)?;
    Array2::from_shape_vec((result.len(), 3), result.into_iter().flatten().collect())
        .map(|array| array.into_pyarray(py).unbind())
        .map_err(value_error)
}

fn triplets(values: PyReadonlyArray2<'_, f64>, shape_error: &str) -> PyResult<Vec<[f64; 3]>> {
    if values.shape().get(1) != Some(&3) {
        return Err(PyValueError::new_err(shape_error.to_owned()));
    }
    values
        .as_slice()
        .map_err(|_| PyValueError::new_err(CONTIGUOUS))?;
    Ok(values
        .as_array()
        .outer_iter()
        .map(|row| [row[0], row[1], row[2]])
        .collect())
}

fn periodic_points(
    values: PyReadonlyArray2<'_, f64>,
) -> PyResult<Vec<Vec<molframe::PeriodicAngle>>> {
    values
        .as_slice()
        .map_err(|_| PyValueError::new_err(CONTIGUOUS))?;
    values
        .as_array()
        .outer_iter()
        .map(|row| {
            row.iter()
                .copied()
                .map(molframe::PeriodicAngle::from_radians)
                .collect::<Result<Vec<_>, _>>()
                .map_err(value_error)
        })
        .collect()
}

fn value_error(error: impl std::fmt::Debug) -> PyErr {
    PyValueError::new_err(format!("{error:?}"))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyMeanSquaredDisplacement>()?;
    module.add_class::<PyTrajectoryDielectricOptions>()?;
    module.add_class::<PyDielectricEstimate>()?;
    module.add_function(wrap_pyfunction!(kmeans, module)?)?;
    module.add_function(wrap_pyfunction!(mean_squared_displacement, module)?)?;
    module.add_function(wrap_pyfunction!(
        trajectory_dielectric_from_dipoles,
        module
    )?)?;
    module.add_function(wrap_pyfunction!(rmsd_to_reference, module)?)?;
    module.add_function(wrap_pyfunction!(pairwise_fitted_rmsd, module)?)?;
    module.add_function(wrap_pyfunction!(pairwise_torus_distance, module)?)?;
    module.add_function(wrap_pyfunction!(generalized_procrustes_mean, module)?)?;
    Ok(())
}
