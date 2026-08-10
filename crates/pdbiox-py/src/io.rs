//! Explicit native read policy and successful diagnostic preservation.

use crate::errors::read_error;
use crate::structure::PyStructure;
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use std::collections::BTreeMap;
use std::path::PathBuf;

#[pyclass(name = "Format", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyFormat {
    Auto,
    Mmcif,
    Pdbml,
    BinaryCif,
    Mmtf,
    Pdb,
    Pqr,
    Pdbqt,
}

#[pyclass(name = "ParseMode", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyParseMode {
    Strict,
    Permissive,
    Recover,
}

#[pyclass(name = "MissingElementPolicy", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyMissingElementPolicy {
    PreserveUnknown,
    InferFromAtomName,
}

#[pyclass(
    name = "AmbiguousResidueBoundaryPolicy",
    frozen,
    eq,
    eq_int,
    from_py_object
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyAmbiguousResidueBoundaryPolicy {
    Reject,
    InferFromFileOrder,
}

#[pyclass(name = "Limits", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyLimits(pdbiox::Limits);

#[pyclass(name = "ReadScope", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyReadScope {
    only_first_model: bool,
    only_atomic_coords: bool,
    discard_hydrogens: bool,
}

#[pymethods]
impl PyReadScope {
    #[new]
    fn new(only_first_model: bool, only_atomic_coords: bool, discard_hydrogens: bool) -> Self {
        Self {
            only_first_model,
            only_atomic_coords,
            discard_hydrogens,
        }
    }

    #[staticmethod]
    fn all() -> Self {
        Self::new(false, false, false)
    }
}

#[pymethods]
impl PyLimits {
    #[new]
    fn new(
        decompressed_bytes: u64,
        compression_ratio: u64,
        rows_per_category: u64,
        nesting_depth: u32,
        dictionary_entries: u32,
    ) -> Self {
        Self(pdbiox::Limits {
            decompressed_bytes,
            compression_ratio,
            rows_per_category,
            nesting_depth,
            dictionary_entries,
        })
    }

    #[staticmethod]
    fn standard() -> Self {
        Self(pdbiox::Limits::default())
    }
}

#[pyclass(name = "ReadOptions", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyReadOptions(pub(crate) pdbiox::ReadOptions);

#[pymethods]
impl PyReadOptions {
    #[new]
    fn new(
        format: PyFormat,
        mode: PyParseMode,
        scope: PyReadScope,
        missing_element_policy: PyMissingElementPolicy,
        ambiguous_residue_boundary_policy: PyAmbiguousResidueBoundaryPolicy,
        limits: PyLimits,
    ) -> Self {
        Self(
            pdbiox::ReadOptions::new()
                .format(format.into())
                .mode(mode.into())
                .only_first_model(scope.only_first_model)
                .only_atomic_coords(scope.only_atomic_coords)
                .discard_hydrogens(scope.discard_hydrogens)
                .missing_element_policy(missing_element_policy.into())
                .ambiguous_residue_boundary_policy(ambiguous_residue_boundary_policy.into())
                .limits(limits.0),
        )
    }

    #[staticmethod]
    fn standard() -> Self {
        Self(pdbiox::ReadOptions::new())
    }
}

#[pyclass(name = "ReadReport", frozen, skip_from_py_object)]
pub(crate) struct PyReadReport {
    structure: PyStructure,
    #[pyo3(get)]
    findings: Vec<String>,
}

#[pyclass(name = "PdbWriteOptions", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPdbWriteOptions(pdbiox::PdbOptions);

#[pymethods]
impl PyPdbWriteOptions {
    #[new]
    fn new(chain_map: BTreeMap<String, String>, hybrid36: bool) -> Self {
        let mut options = pdbiox::PdbOptions::new().hybrid36(hybrid36);
        for (source, target) in chain_map {
            options = options.chain_map(source, target);
        }
        Self(options)
    }
}

#[pymethods]
impl PyReadReport {
    #[getter]
    fn structure(&self) -> PyStructure {
        self.structure.clone()
    }
}

#[pyfunction]
pub(crate) fn read_with_options(
    py: Python<'_>,
    path: PathBuf,
    options: &PyReadOptions,
) -> PyResult<PyReadReport> {
    let options = options.0.clone();
    py.detach(move || pdbiox::read_with_options(path, &options))
        .map(|(structure, findings)| PyReadReport {
            structure: PyStructure::new(structure),
            findings: findings
                .into_iter()
                .map(|value| value.to_string())
                .collect(),
        })
        .map_err(|findings| read_error(py, &findings))
}

#[pyfunction]
#[pyo3(signature = (data, options, *, name=None))]
pub(crate) fn read_bytes(
    py: Python<'_>,
    data: &Bound<'_, PyBytes>,
    options: &PyReadOptions,
    name: Option<String>,
) -> PyResult<PyReadReport> {
    let data = data.as_bytes().to_vec();
    let options = options.0.clone();
    py.detach(move || pdbiox::read_bytes(data, name.as_deref(), &options))
        .map(|(structure, findings)| PyReadReport {
            structure: PyStructure::new(structure),
            findings: findings
                .into_iter()
                .map(|value| value.to_string())
                .collect(),
        })
        .map_err(|findings| read_error(py, &findings))
}

#[pyfunction(signature = (
    structure,
    *,
    block_id = None,
    generate_connection_ids = false,
    connection_type = None
))]
pub(crate) fn write_mmcif(
    py: Python<'_>,
    structure: &PyStructure,
    block_id: Option<String>,
    generate_connection_ids: bool,
    connection_type: Option<String>,
) -> PyResult<String> {
    let structure = structure.structure().clone();
    let mut options = pdbiox::CifWriteOptions::new();
    if let Some(block_id) = block_id {
        options = options.with_block_id(block_id);
    }
    if generate_connection_ids {
        options = options.with_generated_connection_ids();
    }
    if let Some(connection_type) = connection_type {
        options = options.with_connection_type_id(connection_type);
    }
    py.detach(move || pdbiox::write_mmcif_with_options(&structure, &options))
        .map_err(|error| crate::errors::cif_write_error(&error))
}

#[pyfunction]
pub(crate) fn write_bcif<'py>(
    py: Python<'py>,
    structure: &PyStructure,
) -> PyResult<Bound<'py, PyBytes>> {
    let structure = structure.structure().clone();
    py.detach(move || pdbiox::write_bcif(&structure))
        .map(|bytes| PyBytes::new(py, &bytes))
        .map_err(|findings| read_error(py, &findings))
}

#[pyfunction]
pub(crate) fn write_pdb(
    py: Python<'_>,
    structure: &PyStructure,
    options: &PyPdbWriteOptions,
) -> PyResult<String> {
    let structure = structure.structure().clone();
    let options = options.0.clone();
    py.detach(move || pdbiox::write_pdb(&structure, &options))
        .map_err(|findings| read_error(py, &findings))
}

impl From<PyFormat> for pdbiox::Format {
    fn from(value: PyFormat) -> Self {
        match value {
            PyFormat::Auto => Self::Auto,
            PyFormat::Mmcif => Self::Mmcif,
            PyFormat::Pdbml => Self::Pdbml,
            PyFormat::BinaryCif => Self::BinaryCif,
            PyFormat::Mmtf => Self::Mmtf,
            PyFormat::Pdb => Self::Pdb,
            PyFormat::Pqr => Self::Pqr,
            PyFormat::Pdbqt => Self::Pdbqt,
        }
    }
}

impl From<PyParseMode> for pdbiox::ParseMode {
    fn from(value: PyParseMode) -> Self {
        match value {
            PyParseMode::Strict => Self::Strict,
            PyParseMode::Permissive => Self::Permissive,
            PyParseMode::Recover => Self::Recover,
        }
    }
}

impl From<PyMissingElementPolicy> for pdbiox::MissingElementPolicy {
    fn from(value: PyMissingElementPolicy) -> Self {
        match value {
            PyMissingElementPolicy::PreserveUnknown => Self::PreserveUnknown,
            PyMissingElementPolicy::InferFromAtomName => Self::InferFromAtomName,
        }
    }
}

impl From<PyAmbiguousResidueBoundaryPolicy> for pdbiox::AmbiguousResidueBoundaryPolicy {
    fn from(value: PyAmbiguousResidueBoundaryPolicy) -> Self {
        match value {
            PyAmbiguousResidueBoundaryPolicy::Reject => Self::Reject,
            PyAmbiguousResidueBoundaryPolicy::InferFromFileOrder => Self::InferFromFileOrder,
        }
    }
}
