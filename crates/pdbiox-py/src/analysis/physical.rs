//! Explicit polymer and pore geometry without inferred chemistry.

use crate::geometry::coordinates;
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
    let path = coordinates(path)?;
    py.detach(|| pdbiox::analysis::polymer_statistics(&path))
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
pub(crate) struct PyPoreOptions(pdbiox::analysis::PoreProfileOptions);

#[pymethods]
impl PyPoreOptions {
    #[new]
    fn new(
        axis: ([f32; 3], [f32; 3]),
        start: f32,
        end: f32,
        samples: usize,
        search_radius: f32,
        grid_spacing: f32,
        probe_radius: f32,
    ) -> Self {
        Self(pdbiox::analysis::PoreProfileOptions {
            axis_origin: axis.0,
            axis_direction: axis.1,
            start,
            end,
            samples,
            search_radius,
            grid_spacing,
            probe_radius,
        })
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

#[pyfunction]
pub(crate) fn pore_profile(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    radii: PyReadonlyArray1<'_, f32>,
    options: &PyPoreOptions,
) -> PyResult<Vec<PyPoreSample>> {
    let positions = coordinates(positions)?;
    let radius_values = radii.as_slice()?.to_vec();
    drop(radii);
    let options = options.0;
    py.detach(|| pdbiox::analysis::pore_profile(&positions, &radius_values, options))
        .map(|samples| {
            samples
                .into_iter()
                .map(|value| PyPoreSample {
                    axial_coordinate: value.axial_coordinate,
                    centre: value.centre,
                    radius: value.radius,
                })
                .collect()
        })
        .map_err(value_error)
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
