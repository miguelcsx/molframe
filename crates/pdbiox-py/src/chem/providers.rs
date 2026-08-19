//! Python handles for public chemistry providers and free reference queries.

use super::components::PyComponent;
use super::{PyComponentDictionary, PyElement, PyElementProperties, PyIonicRadius, PyRadiusSet};
use ::pdbiox;
use pyo3::prelude::*;
use pyo3::types::PyList;
use std::path::PathBuf;
use std::sync::Arc;

#[pyclass(name = "MemoryProvider", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyMemoryProvider(pub(crate) Arc<pdbiox::MemoryProvider>);

/// A cheap owned provider handle used by free facade operations.
///
/// The enum is resolved once at the Python boundary; all structure-wide work
/// still receives the native `ComponentProvider` trait object and runs in Rust.
#[derive(Clone)]
pub(crate) enum PyComponentProvider {
    Cif(Arc<pdbiox::CifProvider>),
    Memory(Arc<pdbiox::MemoryProvider>),
}

impl pdbiox::ComponentProvider for PyComponentProvider {
    fn get(
        &self,
        component_id: &str,
    ) -> Result<Option<Arc<pdbiox::Component>>, pdbiox::Diagnostic> {
        match self {
            Self::Cif(provider) => provider.get(component_id),
            Self::Memory(provider) => provider.get(component_id),
        }
    }

    fn version(&self) -> &pdbiox::DictionaryVersion {
        match self {
            Self::Cif(provider) => provider.version(),
            Self::Memory(provider) => provider.version(),
        }
    }
}

pub(crate) fn extract_provider(provider: &Bound<'_, PyAny>) -> PyResult<PyComponentProvider> {
    if let Ok(value) = provider.extract::<PyRef<'_, PyComponentDictionary>>() {
        return Ok(PyComponentProvider::Cif(value.0.clone()));
    }
    if let Ok(value) = provider.extract::<PyRef<'_, PyMemoryProvider>>() {
        return Ok(PyComponentProvider::Memory(value.0.clone()));
    }
    Err(pyo3::exceptions::PyTypeError::new_err(
        "provider must be ComponentDictionary/CifProvider or MemoryProvider",
    ))
}

#[pymethods]
impl PyMemoryProvider {
    #[new]
    fn new(version: &str, components: &Bound<'_, PyList>) -> PyResult<Self> {
        let components = components
            .iter()
            .map(|value| {
                value
                    .extract::<PyRef<'_, PyComponent>>()
                    .map_err(PyErr::from)
                    .map(|component| (*component.0).clone())
            })
            .collect::<PyResult<Vec<_>>>()?;
        pdbiox::MemoryProvider::new(pdbiox::DictionaryVersion::new(version), components)
            .map(|provider| Self(Arc::new(provider)))
            .map_err(value_error)
    }

    #[getter]
    fn version(&self) -> &str {
        pdbiox::ComponentProvider::version(self.0.as_ref()).as_str()
    }

    fn get(&self, component_id: &str) -> PyResult<Option<PyComponent>> {
        pdbiox::ComponentProvider::get(self.0.as_ref(), component_id)
            .map(|component| component.map(PyComponent))
            .map_err(value_error)
    }
}

#[pyfunction]
pub(crate) fn read_ccd(
    py: Python<'_>,
    path: PathBuf,
    version: &str,
) -> PyResult<(PyComponentDictionary, Vec<String>)> {
    let bytes = std::fs::read(&path).map_err(|error| {
        pyo3::exceptions::PyOSError::new_err(format!("cannot read CCD {}: {error}", path.display()))
    })?;
    let version = pdbiox::DictionaryVersion::new(version);
    py.detach(move || {
        let input = pdbiox::InputBuffer::from_bytes(bytes);
        pdbiox::read_ccd(&input, version)
    })
    .map(|(provider, findings)| {
        (
            PyComponentDictionary(Arc::new(provider)),
            findings
                .into_iter()
                .map(|finding| finding.to_string())
                .collect(),
        )
    })
    .map_err(|findings| crate::errors::read_error(py, &findings))
}

#[pyfunction]
pub(crate) fn element_properties(element: &PyElement) -> Option<PyElementProperties> {
    pdbiox::element_properties(element.0).map(PyElementProperties::from)
}

#[pyfunction]
pub(crate) fn vdw_radius(element: &PyElement, set: PyRadiusSet) -> Option<f32> {
    pdbiox::vdw_radius(element.0, set.into())
}

#[pyfunction]
pub(crate) fn ionic_radii(element: &PyElement) -> PyResult<Vec<PyIonicRadius>> {
    pdbiox::ionic_radii(element.0)
        .map(|values| values.iter().cloned().map(PyIonicRadius::from).collect())
        .map_err(value_error)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyMemoryProvider>()?;
    module.add_function(wrap_pyfunction!(read_ccd, module)?)?;
    module.add_function(wrap_pyfunction!(element_properties, module)?)?;
    module.add_function(wrap_pyfunction!(vdw_radius, module)?)?;
    module.add_function(wrap_pyfunction!(ionic_radii, module)?)?;
    Ok(())
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    pyo3::exceptions::PyValueError::new_err(error.to_string())
}
