//! Explicit native read policy and successful diagnostic preservation.

use crate::contract::PyDiagnostic;
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
pub(crate) struct PyLimits(pub(crate) molframe::Limits);

impl PyLimits {
    pub(crate) const fn from_inner(value: molframe::Limits) -> Self {
        Self(value)
    }
}

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
        Self(molframe::Limits {
            decompressed_bytes,
            compression_ratio,
            rows_per_category,
            nesting_depth,
            dictionary_entries,
        })
    }

    #[staticmethod]
    fn standard() -> Self {
        Self(molframe::Limits::default())
    }
}

#[pyclass(name = "ReadOptions", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyReadOptions(pub(crate) molframe::ReadOptions);

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
            molframe::ReadOptions::new()
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
        Self(molframe::ReadOptions::new())
    }
}

#[pyclass(name = "ReadReport", frozen, skip_from_py_object)]
pub(crate) struct PyReadReport {
    pub(crate) structure: PyStructure,
    #[pyo3(get)]
    pub(crate) findings: Vec<PyDiagnostic>,
}

#[pyclass(name = "PdbWriteOptions", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPdbWriteOptions(pub(crate) molframe::PdbOptions);

/// Identifier namespace projected into fixed-width PDB fields.
#[pyclass(name = "PdbIdentifierNamespace", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum PyPdbIdentifierNamespace {
    #[default]
    Label,
    Auth,
}

impl From<PyPdbIdentifierNamespace> for molframe::PdbIdentifierNamespace {
    fn from(value: PyPdbIdentifierNamespace) -> Self {
        match value {
            PyPdbIdentifierNamespace::Label => Self::Label,
            PyPdbIdentifierNamespace::Auth => Self::Auth,
        }
    }
}

#[pymethods]
impl PyPdbWriteOptions {
    #[new]
    #[pyo3(signature = (chain_map=None, hybrid36=false, namespace=PyPdbIdentifierNamespace::Label))]
    fn new(
        chain_map: Option<BTreeMap<String, String>>,
        hybrid36: bool,
        namespace: PyPdbIdentifierNamespace,
    ) -> Self {
        let mut options = molframe::PdbOptions::new()
            .hybrid36(hybrid36)
            .namespace(namespace.into());
        for (source, target) in chain_map.into_iter().flatten() {
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
    py.detach(move || molframe::read_with_options(path, &options))
        .map(|(structure, findings)| PyReadReport {
            structure: PyStructure::new(structure),
            findings: findings.into_iter().map(Into::into).collect(),
        })
        .map_err(|findings| read_error(py, &findings))
}

#[pyfunction]
pub(crate) fn read_with_diagnostics(py: Python<'_>, path: PathBuf) -> PyResult<PyReadReport> {
    py.detach(move || molframe::read_with_diagnostics(path))
        .map(|(structure, findings)| PyReadReport {
            structure: PyStructure::new(structure),
            findings: findings.into_iter().map(Into::into).collect(),
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
    py.detach(move || molframe::read_bytes(data, name.as_deref(), &options))
        .map(|(structure, findings)| PyReadReport {
            structure: PyStructure::new(structure),
            findings: findings.into_iter().map(Into::into).collect(),
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
    let mut options = molframe::CifWriteOptions::new();
    if let Some(block_id) = block_id {
        options = options.with_block_id(block_id);
    }
    if generate_connection_ids {
        options = options.with_generated_connection_ids();
    }
    if let Some(connection_type) = connection_type {
        options = options.with_connection_type_id(connection_type);
    }
    py.detach(move || molframe::write_mmcif_with_options(&structure, &options))
        .map_err(|error| crate::errors::cif_write_error(&error))
}

#[pyfunction]
pub(crate) fn write_bcif<'py>(
    py: Python<'py>,
    structure: &PyStructure,
) -> PyResult<Bound<'py, PyBytes>> {
    let structure = structure.structure().clone();
    py.detach(move || molframe::write_bcif(&structure))
        .map(|bytes| PyBytes::new(py, &bytes))
        .map_err(|findings| read_error(py, &findings))
}

#[pyfunction(signature = (
    structure,
    *,
    block_id = None,
    generate_connection_ids = false,
    connection_type = None
))]
pub(crate) fn write_bcif_with_options<'py>(
    py: Python<'py>,
    structure: &PyStructure,
    block_id: Option<String>,
    generate_connection_ids: bool,
    connection_type: Option<String>,
) -> PyResult<Bound<'py, PyBytes>> {
    let structure = structure.structure().clone();
    let mut options = molframe::CifWriteOptions::new();
    if let Some(block_id) = block_id {
        options = options.with_block_id(block_id);
    }
    if generate_connection_ids {
        options = options.with_generated_connection_ids();
    }
    if let Some(connection_type) = connection_type {
        options = options.with_connection_type_id(connection_type);
    }
    py.detach(move || molframe::write_bcif_with_options(&structure, &options))
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
    py.detach(move || molframe::write_pdb(&structure, &options))
        .map_err(|findings| read_error(py, &findings))
}

#[pyfunction]
pub(crate) fn write(py: Python<'_>, path: PathBuf, structure: &PyStructure) -> PyResult<()> {
    let structure = structure.structure().clone();
    py.detach(move || molframe::write(path, &structure))
        .map_err(|findings| read_error(py, &findings))
}

#[pyfunction]
#[pyo3(signature = (data, options=None))]
pub(crate) fn read_mmtf(
    py: Python<'_>,
    data: &Bound<'_, PyBytes>,
    options: Option<&PyReadOptions>,
) -> PyResult<PyReadReport> {
    read_variant(py, data, options, molframe::pdb::read_mmtf)
}

#[pyfunction]
#[pyo3(signature = (data, options=None))]
pub(crate) fn read_pdb(
    py: Python<'_>,
    data: &Bound<'_, PyAny>,
    options: Option<&PyReadOptions>,
) -> PyResult<PyReadReport> {
    let data = match data.cast::<PyBytes>() {
        Ok(bytes) => bytes.as_bytes().to_vec(),
        Err(_) => data.extract::<String>()?.into_bytes(),
    };
    let options = options.map_or_else(molframe::ReadOptions::new, |value| value.0.clone());
    py.detach(move || {
        let input = molframe::InputBuffer::from_bytes(data);
        molframe::pdb::read(&input, &options)
    })
    .map(|(structure, findings)| PyReadReport {
        structure: PyStructure::new(structure),
        findings: findings.into_iter().map(Into::into).collect(),
    })
    .map_err(|findings| read_error(py, &findings))
}

#[pyfunction]
#[pyo3(signature = (data, options=None))]
pub(crate) fn read_pqr(
    py: Python<'_>,
    data: &Bound<'_, PyBytes>,
    options: Option<&PyReadOptions>,
) -> PyResult<PyReadReport> {
    read_variant(py, data, options, molframe::pdb::read_pqr)
}

#[pyfunction]
#[pyo3(signature = (data, options=None))]
pub(crate) fn read_pdbqt(
    py: Python<'_>,
    data: &Bound<'_, PyBytes>,
    options: Option<&PyReadOptions>,
) -> PyResult<PyReadReport> {
    read_variant(py, data, options, molframe::pdb::read_pdbqt)
}

fn read_variant(
    py: Python<'_>,
    data: &Bound<'_, PyBytes>,
    options: Option<&PyReadOptions>,
    reader: fn(&molframe::InputBuffer, &molframe::ReadOptions) -> molframe::ReadResult,
) -> PyResult<PyReadReport> {
    let data = data.as_bytes().to_vec();
    let options = options.map_or_else(molframe::ReadOptions::new, |value| value.0.clone());
    py.detach(move || {
        let input = molframe::InputBuffer::from_bytes(data);
        reader(&input, &options)
    })
    .map(|(structure, findings)| PyReadReport {
        structure: PyStructure::new(structure),
        findings: findings.into_iter().map(Into::into).collect(),
    })
    .map_err(|findings| read_error(py, &findings))
}

#[pyfunction]
pub(crate) fn write_mmtf<'py>(
    py: Python<'py>,
    structure: &PyStructure,
) -> PyResult<Bound<'py, PyBytes>> {
    let structure = structure.structure().clone();
    py.detach(move || molframe::write_mmtf(&structure))
        .map(|bytes| PyBytes::new(py, &bytes))
        .map_err(|findings| read_error(py, &findings))
}

#[pyfunction]
pub(crate) fn write_pqr(
    py: Python<'_>,
    structure: &PyStructure,
    options: &PyPdbWriteOptions,
) -> PyResult<String> {
    let structure = structure.structure().clone();
    let options = options.0.clone();
    py.detach(move || molframe::write_pqr(&structure, &options))
        .map_err(|findings| read_error(py, &findings))
}

#[pyfunction]
pub(crate) fn write_pdbqt(
    py: Python<'_>,
    structure: &PyStructure,
    options: &PyPdbWriteOptions,
) -> PyResult<String> {
    let structure = structure.structure().clone();
    let options = options.0.clone();
    py.detach(move || molframe::write_pdbqt(&structure, &options))
        .map_err(|findings| read_error(py, &findings))
}

impl From<PyFormat> for molframe::Format {
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

impl From<PyParseMode> for molframe::ParseMode {
    fn from(value: PyParseMode) -> Self {
        match value {
            PyParseMode::Strict => Self::Strict,
            PyParseMode::Permissive => Self::Permissive,
            PyParseMode::Recover => Self::Recover,
        }
    }
}

impl From<PyMissingElementPolicy> for molframe::MissingElementPolicy {
    fn from(value: PyMissingElementPolicy) -> Self {
        match value {
            PyMissingElementPolicy::PreserveUnknown => Self::PreserveUnknown,
            PyMissingElementPolicy::InferFromAtomName => Self::InferFromAtomName,
        }
    }
}

impl From<PyAmbiguousResidueBoundaryPolicy> for molframe::AmbiguousResidueBoundaryPolicy {
    fn from(value: PyAmbiguousResidueBoundaryPolicy) -> Self {
        match value {
            PyAmbiguousResidueBoundaryPolicy::Reject => Self::Reject,
            PyAmbiguousResidueBoundaryPolicy::InferFromFileOrder => Self::InferFromFileOrder,
        }
    }
}
