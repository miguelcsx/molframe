//! Python connected-component operations for indexed surfaces.

use super::types::PyIndexedSurfaceMesh;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

#[pyclass(name = "SurfaceComponent", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySurfaceComponent {
    #[pyo3(get)]
    faces: Vec<u32>,
    #[pyo3(get)]
    area: f64,
}

impl From<pdbiox::surface::SurfaceComponent> for PySurfaceComponent {
    fn from(value: pdbiox::surface::SurfaceComponent) -> Self {
        Self {
            faces: value.faces,
            area: value.area,
        }
    }
}

#[pyclass(name = "SurfaceComponentFilter", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PySurfaceComponentFilter {
    native: pdbiox::surface::SurfaceComponentFilter,
}

#[pymethods]
impl PySurfaceComponentFilter {
    #[new]
    #[pyo3(signature = (minimum_area=0.0, maximum_components=None))]
    fn new(minimum_area: f64, maximum_components: Option<usize>) -> Self {
        Self {
            native: pdbiox::surface::SurfaceComponentFilter {
                minimum_area,
                maximum_components,
            },
        }
    }

    #[getter]
    fn minimum_area(&self) -> f64 {
        self.native.minimum_area
    }
    #[getter]
    fn maximum_components(&self) -> Option<usize> {
        self.native.maximum_components
    }
}

#[pyclass(name = "SurfaceComponentError", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PySurfaceComponentError {
    InvalidFilter,
    MeshTooLarge,
}

#[pyfunction]
fn surface_components(mesh: &PyIndexedSurfaceMesh) -> Vec<PySurfaceComponent> {
    pdbiox::surface::surface_components(&mesh.native)
        .into_iter()
        .map(Into::into)
        .collect()
}

#[pyfunction]
fn filter_surface_components(
    mesh: &PyIndexedSurfaceMesh,
    filter: &PySurfaceComponentFilter,
) -> PyResult<PyIndexedSurfaceMesh> {
    pdbiox::surface::filter_surface_components(&mesh.native, filter.native)
        .map(PyIndexedSurfaceMesh::from_native)
        .map_err(|error| PyValueError::new_err(error.to_string()))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PySurfaceComponent>()?;
    module.add_class::<PySurfaceComponentFilter>()?;
    module.add_class::<PySurfaceComponentError>()?;
    module.add_function(wrap_pyfunction!(surface_components, module)?)?;
    module.add_function(wrap_pyfunction!(filter_surface_components, module)?)
}
