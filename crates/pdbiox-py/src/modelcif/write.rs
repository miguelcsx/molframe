//! Canonical `ModelCIF` output through the native writer.

use crate::cif_write::{PyCifWriteOptions, write_error};
use crate::modelcif::PyModelCif;
use crate::structure::PyStructure;
use pyo3::prelude::*;

#[pyfunction(name = "modelcif_write_canonical")]
#[pyo3(signature = (structure, model_cif, options=None))]
pub(crate) fn write_canonical(
    py: Python<'_>,
    structure: &PyStructure,
    model_cif: &PyModelCif,
    options: Option<&PyCifWriteOptions>,
) -> PyResult<String> {
    write_native(py, structure, model_cif, options)
}

#[pyfunction(name = "modelcif_write_canonical_with_options")]
pub(crate) fn write_canonical_with_options(
    py: Python<'_>,
    structure: &PyStructure,
    model_cif: &PyModelCif,
    options: &PyCifWriteOptions,
) -> PyResult<String> {
    write_native(py, structure, model_cif, Some(options))
}

#[pyfunction(name = "modelcif_attach")]
pub(crate) fn attach(
    py: Python<'_>,
    structure: &PyStructure,
    model_cif: &PyModelCif,
) -> PyResult<PyStructure> {
    let structure = structure.structure().clone();
    let model_cif = model_cif.inner().clone();
    py.detach(move || {
        Ok(PyStructure::new(
            structure.with_extension(pdbiox::MODEL_CIF_EXTENSION, model_cif),
        ))
    })
}

fn write_native(
    py: Python<'_>,
    structure: &PyStructure,
    model_cif: &PyModelCif,
    options: Option<&PyCifWriteOptions>,
) -> PyResult<String> {
    let structure = structure.structure().clone();
    let model_cif = std::sync::Arc::clone(&model_cif.inner);
    let options = options.map_or_else(PyCifWriteOptions::standard, Clone::clone);
    py.detach(move || {
        pdbiox::modelcif::write_canonical_with_options(&structure, &model_cif, &options.native())
    })
    .map_err(write_error)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(write_canonical, module)?)?;
    module.add_function(wrap_pyfunction!(write_canonical_with_options, module)?)?;
    module.add_function(wrap_pyfunction!(attach, module)?)?;
    module.add("MODEL_CIF_EXTENSION", pdbiox::MODEL_CIF_EXTENSION)?;
    Ok(())
}
