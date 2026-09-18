//! Explicit polymer and pore geometry without inferred chemistry.

use crate::geometry::borrowed_coordinates;
use numpy::{PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "PolymerStatistics", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyPolymerStatistics {
    #[pyo3(get)]
    contour_length: f64,
    #[pyo3(get)]
    end_to_end_distance: f64,
    #[pyo3(get)]
    mean_segment_length: f64,
    #[pyo3(get)]
    adjacent_tangent_correlation: Option<f64>,
    #[pyo3(get)]
    persistence_length: Option<f64>,
}

#[pyfunction]
pub(crate) fn polymer_statistics(
    py: Python<'_>,
    path: PyReadonlyArray2<'_, f32>,
) -> PyResult<PyPolymerStatistics> {
    let path = borrowed_coordinates(&path)?;
    py.detach(|| molframe::analysis::polymer_statistics(path))
        .map(|value| PyPolymerStatistics {
            contour_length: value.contour_length,
            end_to_end_distance: value.end_to_end_distance,
            mean_segment_length: value.mean_segment_length,
            adjacent_tangent_correlation: value.adjacent_tangent_correlation,
            persistence_length: value.persistence_length,
        })
        .map_err(value_error)
}

#[pyclass(name = "PoreOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyPoreOptions(molframe::analysis::PoreProfileOptions);

#[pymethods]
impl PyPoreOptions {
    #[new]
    #[pyo3(signature = (axis, start, end, samples, search_radius, grid_spacing, probe_radius, *, memory_limit_bytes = 100_000_000))]
    fn new(
        axis: ([f32; 3], [f32; 3]),
        start: f32,
        end: f32,
        samples: usize,
        search_radius: f32,
        grid_spacing: f32,
        probe_radius: f32,
        memory_limit_bytes: usize,
    ) -> Self {
        Self(molframe::analysis::PoreProfileOptions {
            axis_origin: axis.0,
            axis_direction: axis.1,
            start,
            end,
            samples,
            search_radius,
            grid_spacing,
            probe_radius,
            memory_limit_bytes,
        })
    }
}

impl PyPoreOptions {
    pub(crate) const fn native(&self) -> molframe::analysis::PoreProfileOptions {
        self.0
    }

    pub(crate) const fn from_native(options: molframe::analysis::PoreProfileOptions) -> Self {
        Self(options)
    }
}

#[pyclass(name = "PoreSample", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyPoreSample {
    #[pyo3(get)]
    axial_coordinate: f32,
    #[pyo3(get)]
    centre: [f32; 3],
    #[pyo3(get)]
    radius: f32,
}

impl From<molframe::analysis::PoreSample> for PyPoreSample {
    fn from(value: molframe::analysis::PoreSample) -> Self {
        Self {
            axial_coordinate: value.axial_coordinate,
            centre: value.centre,
            radius: value.radius,
        }
    }
}

pub(crate) fn pore_analysis(
    py: Python<'_>,
    analysis: molframe::Analysis<Vec<molframe::analysis::PoreSample>>,
) -> PyResult<crate::contract::PyAnalysis> {
    crate::contract::analysis_with_value(py, analysis, |py, values| {
        let values = values
            .into_iter()
            .map(|value| Py::new(py, PyPoreSample::from(value)))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(pyo3::types::PyList::new(py, values)?.unbind().into_any())
    })
}

#[pyfunction]
pub(crate) fn pore_profile(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    radii: PyReadonlyArray1<'_, f32>,
    options: &PyPoreOptions,
) -> PyResult<Vec<PyPoreSample>> {
    let positions = borrowed_coordinates(&positions)?;
    let radius_values = radii.as_slice()?;
    let options = options.0;
    py.detach(|| molframe::analysis::pore_profile(positions, radius_values, options))
        .map(|samples| samples.into_iter().map(PyPoreSample::from).collect())
        .map_err(value_error)
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
