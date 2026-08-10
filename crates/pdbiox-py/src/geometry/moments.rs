//! Vectorized shape moments and best-fit planes delegated to Rust.

use super::coordinates;
use numpy::{PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "Axes", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyAxes {
    #[pyo3(get)]
    values: [f64; 3],
    #[pyo3(get)]
    vectors: [[f64; 3]; 3],
}

#[pyclass(name = "Plane", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyPlane {
    #[pyo3(get)]
    centre: [f64; 3],
    #[pyo3(get)]
    normal: [f64; 3],
}

#[pyclass(name = "EigenOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyEigenOptions {
    pub(crate) inner: pdbiox::EigenOptions,
}

#[pymethods]
impl PyEigenOptions {
    #[new]
    const fn new(relative_tolerance: f64, maximum_sweeps: usize) -> Self {
        Self {
            inner: pdbiox::EigenOptions {
                relative_tolerance,
                maximum_sweeps,
            },
        }
    }

    #[getter]
    const fn relative_tolerance(&self) -> f64 {
        self.inner.relative_tolerance
    }

    #[getter]
    const fn maximum_sweeps(&self) -> usize {
        self.inner.maximum_sweeps
    }
}

#[pyfunction]
pub(crate) fn centroid(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
) -> PyResult<Option<[f64; 3]>> {
    let positions = coordinates(positions)?;
    Ok(py.detach(|| pdbiox::centroid(&positions)))
}

#[pyfunction]
pub(crate) fn centre_of_mass(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    masses: &Bound<'_, PyAny>,
) -> PyResult<Option<[f64; 3]>> {
    let positions = coordinates(positions)?;
    let mass_values = mass_values(masses, positions.len())?;
    Ok(py.detach(|| pdbiox::centre_of_mass(&positions, &mass_values)))
}

#[pyfunction]
pub(crate) fn radius_of_gyration(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    masses: &Bound<'_, PyAny>,
) -> PyResult<Option<f64>> {
    let positions = coordinates(positions)?;
    let mass_values = mass_values(masses, positions.len())?;
    Ok(py.detach(|| pdbiox::radius_of_gyration(&positions, &mass_values)))
}

#[pyfunction]
pub(crate) fn inertia_tensor(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    masses: &Bound<'_, PyAny>,
) -> PyResult<Option<[[f64; 3]; 3]>> {
    let positions = coordinates(positions)?;
    let mass_values = mass_values(masses, positions.len())?;
    Ok(py.detach(|| pdbiox::inertia_tensor(&positions, &mass_values)))
}

#[pyfunction]
pub(crate) fn principal_axes(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    masses: &Bound<'_, PyAny>,
    options: &PyEigenOptions,
) -> PyResult<Option<PyAxes>> {
    let positions = coordinates(positions)?;
    let mass_values = mass_values(masses, positions.len())?;
    py.detach(|| pdbiox::principal_axes_with_options(&positions, &mass_values, options.inner))
        .map(|value| value.map(PyAxes::from))
        .map_err(eigen_error)
}

#[pyfunction]
pub(crate) fn asphericity(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    options: &PyEigenOptions,
) -> PyResult<Option<f64>> {
    let positions = coordinates(positions)?;
    py.detach(|| pdbiox::asphericity_with_options(&positions, options.inner))
        .map_err(eigen_error)
}

#[pyfunction]
pub(crate) fn gyration_axes(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    options: &PyEigenOptions,
) -> PyResult<Option<PyAxes>> {
    let positions = coordinates(positions)?;
    py.detach(|| pdbiox::gyration_axes_with_options(&positions, options.inner))
        .map(|value| value.map(PyAxes::from))
        .map_err(eigen_error)
}

#[pyfunction]
pub(crate) fn best_fit_plane(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    options: &PyEigenOptions,
) -> PyResult<Option<PyPlane>> {
    let positions = coordinates(positions)?;
    py.detach(|| pdbiox::best_fit_plane_with_options(&positions, options.inner))
        .map(|value| value.map(PyPlane::from))
        .map_err(eigen_error)
}

#[pyfunction]
pub(crate) fn plane_deviation(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    options: &PyEigenOptions,
) -> PyResult<Option<f64>> {
    let positions = coordinates(positions)?;
    py.detach(|| pdbiox::plane_deviation_with_options(&positions, options.inner))
        .map_err(eigen_error)
}

impl From<pdbiox::Decomposition<3>> for PyAxes {
    fn from(value: pdbiox::Decomposition<3>) -> Self {
        Self {
            values: value.values,
            vectors: value.vectors,
        }
    }
}

impl From<pdbiox::Plane> for PyPlane {
    fn from(value: pdbiox::Plane) -> Self {
        Self {
            centre: value.centre,
            normal: value.normal,
        }
    }
}

fn mass_values(value: &Bound<'_, PyAny>, positions: usize) -> PyResult<Vec<f64>> {
    if value.is_none() {
        return Ok(Vec::new());
    }
    let array = value.extract::<PyReadonlyArray1<'_, f64>>()?;
    let masses = array.as_slice()?.to_vec();
    if masses.len() != positions {
        return Err(PyValueError::new_err(
            "masses must contain exactly one value per position",
        ));
    }
    Ok(masses)
}

fn eigen_error(error: pdbiox::EigenError) -> PyErr {
    PyValueError::new_err(format!("eigendecomposition failed: {error:?}"))
}
