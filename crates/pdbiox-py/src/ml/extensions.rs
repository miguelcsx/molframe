//! Stable Arrow extension metadata exposed without requiring a Python Arrow import.

use pyo3::prelude::*;
use pyo3::types::PyDict;

#[pyclass(name = "ExportCost", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyExportCost {
    ZeroCopy,
    Decode,
    Copy,
}

#[pymethods]
impl PyExportCost {
    #[getter]
    fn label(slf: PyRef<'_, Self>) -> &'static str {
        let value = *slf;
        drop(slf);
        match value {
            Self::ZeroCopy => "zero-copy",
            Self::Decode => "decode",
            Self::Copy => "copy",
        }
    }
}

#[pyclass(name = "PdbioxExtension", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyPdbioxExtension {
    AtomIndex,
    ResidueIndex,
    ChainIndex,
    EntityIndex,
    Coordinates3f,
    Element,
    SymbolId,
    Altloc,
    Selection,
    Validity,
}

#[pymethods]
impl PyPdbioxExtension {
    #[getter]
    fn name(slf: PyRef<'_, Self>) -> &'static str {
        let value = *slf;
        drop(slf);
        value.inner().name()
    }

    #[getter]
    fn storage_type(slf: PyRef<'_, Self>) -> &'static str {
        let value = *slf;
        drop(slf);
        match value {
            Self::AtomIndex
            | Self::ResidueIndex
            | Self::ChainIndex
            | Self::EntityIndex
            | Self::SymbolId
            | Self::Altloc => "uint32",
            Self::Coordinates3f => "fixed_size_list<float32, 3>",
            Self::Element => "uint8",
            Self::Selection => "binary",
            Self::Validity => "bool",
        }
    }
}

impl PyPdbioxExtension {
    const fn inner(self) -> pdbiox::PdbioxExtension {
        match self {
            Self::AtomIndex => pdbiox::PdbioxExtension::AtomIndex,
            Self::ResidueIndex => pdbiox::PdbioxExtension::ResidueIndex,
            Self::ChainIndex => pdbiox::PdbioxExtension::ChainIndex,
            Self::EntityIndex => pdbiox::PdbioxExtension::EntityIndex,
            Self::Coordinates3f => pdbiox::PdbioxExtension::Coordinates3f,
            Self::Element => pdbiox::PdbioxExtension::Element,
            Self::SymbolId => pdbiox::PdbioxExtension::SymbolId,
            Self::Altloc => pdbiox::PdbioxExtension::Altloc,
            Self::Selection => pdbiox::PdbioxExtension::Selection,
            Self::Validity => pdbiox::PdbioxExtension::Validity,
        }
    }
}

#[pyfunction]
pub(crate) fn extension_name(metadata: &Bound<'_, PyDict>) -> Option<String> {
    metadata
        .get_item("ARROW:extension:name")
        .ok()
        .flatten()
        .and_then(|value| value.extract::<String>().ok())
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyExportCost>()?;
    module.add_class::<PyPdbioxExtension>()?;
    module.add_function(wrap_pyfunction!(extension_name, module)?)?;
    Ok(())
}
