//! Explicit CIF writer options and value rendering.

use crate::cif_document::PyCifValue;
use crate::structure::PyStructure;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

pyo3::create_exception!(_native, CifWriteError, PyValueError);

#[pyclass(name = "CifWriteOptions", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyCifWriteOptions {
    block_id: Option<String>,
    generate_connection_ids: bool,
    connection_type: Option<String>,
}

#[pymethods]
impl PyCifWriteOptions {
    #[new]
    #[pyo3(signature = (*, block_id=None, generate_connection_ids=false, connection_type=None))]
    fn new(
        block_id: Option<String>,
        generate_connection_ids: bool,
        connection_type: Option<String>,
    ) -> Self {
        Self {
            block_id,
            generate_connection_ids,
            connection_type,
        }
    }

    #[staticmethod]
    pub(crate) fn standard() -> Self {
        Self::new(None, false, None)
    }

    #[getter]
    fn block_id(&self) -> Option<String> {
        self.block_id.clone()
    }

    #[getter]
    const fn generate_connection_ids(&self) -> bool {
        self.generate_connection_ids
    }

    #[getter]
    fn connection_type(&self) -> Option<String> {
        self.connection_type.clone()
    }

    fn with_block_id(&self, block_id: String) -> Self {
        let mut value = self.clone();
        value.block_id = Some(block_id);
        value
    }

    fn with_generated_connection_ids(&self) -> Self {
        let mut value = self.clone();
        value.generate_connection_ids = true;
        value
    }

    fn with_connection_type_id(&self, value: String) -> Self {
        let mut options = self.clone();
        options.connection_type = Some(value);
        options
    }
}

pub(crate) fn write_error(error: impl std::fmt::Display) -> PyErr {
    CifWriteError::new_err(error.to_string())
}

impl PyCifWriteOptions {
    pub(crate) fn native(&self) -> pdbiox::CifWriteOptions {
        let mut options = pdbiox::CifWriteOptions::new();
        if let Some(block_id) = &self.block_id {
            options = options.with_block_id(block_id.clone());
        }
        if self.generate_connection_ids {
            options = options.with_generated_connection_ids();
        }
        if let Some(connection_type) = &self.connection_type {
            options = options.with_connection_type_id(connection_type.clone());
        }
        options
    }
}

#[pyfunction(name = "write_canonical")]
#[pyo3(signature = (structure, options=None))]
pub(crate) fn write_canonical(
    py: Python<'_>,
    structure: &PyStructure,
    options: Option<&PyCifWriteOptions>,
) -> PyResult<String> {
    write_native(py, structure, options)
}

#[pyfunction(name = "write_canonical_with_options")]
pub(crate) fn write_canonical_with_options(
    py: Python<'_>,
    structure: &PyStructure,
    options: &PyCifWriteOptions,
) -> PyResult<String> {
    write_native(py, structure, Some(options))
}

#[pyfunction]
pub(crate) fn quote_text(py: Python<'_>, text: &str) -> String {
    py.detach(move || -> String { pdbiox::cif::quote_text(text) })
}

#[pyfunction]
pub(crate) fn render_value(py: Python<'_>, value: &PyCifValue) -> String {
    py.detach(move || -> String { pdbiox::cif::render_value(&value.inner) })
}

fn write_native(
    py: Python<'_>,
    structure: &PyStructure,
    options: Option<&PyCifWriteOptions>,
) -> PyResult<String> {
    let structure = structure.structure().clone();
    let options = options.map_or_else(PyCifWriteOptions::standard, Clone::clone);
    py.detach(move || pdbiox::write_mmcif_with_options(&structure, &options.native()))
        .map_err(|error| CifWriteError::new_err(error.to_string()))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("CifWriteError", module.py().get_type::<CifWriteError>())?;
    module.add_class::<PyCifWriteOptions>()?;
    module.add_function(wrap_pyfunction!(write_canonical, module)?)?;
    module.add_function(wrap_pyfunction!(write_canonical_with_options, module)?)?;
    module.add_function(wrap_pyfunction!(quote_text, module)?)?;
    module.add_function(wrap_pyfunction!(render_value, module)?)?;
    Ok(())
}
