//! Central facade handles and stateless namespace functions.

use crate::hierarchy::{PyAtoms, PyChains, PyModels, PyResidueSelection, PyResidues};
#[cfg(feature = "geometry")]
use numpy::PyReadonlyArray2;
#[cfg(feature = "geometry")]
use numpy::PyUntypedArrayMethods;
use numpy::ndarray::ArrayView2;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods};
use pyo3::prelude::*;
use pyo3::pybacked::PyBackedBytes;
use pyo3::types::PyAny;
use std::fmt;
use std::path::PathBuf;

#[cfg(feature = "analysis")]
mod analysis;
#[cfg(feature = "analysis")]
mod bonds;
#[cfg(feature = "analysis")]
pub(crate) use analysis::{PyContactTable, atom_contacts, contacts};

#[derive(Clone, Debug)]
#[pyclass(name = "Structure", frozen, skip_from_py_object)]
pub(crate) struct PyStructure {
    pub(crate) inner: molframe::Structure,
}
#[derive(Clone, Debug)]
#[pyclass(frozen, skip_from_py_object)]
struct SnapshotOwner {
    structure: molframe::Structure,
}
impl PyStructure {
    pub(crate) const fn new(inner: molframe::Structure) -> Self {
        Self { inner }
    }

    fn coordinates_view<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let owner = Bound::new(
            py,
            SnapshotOwner {
                structure: self.inner.clone(),
            },
        )?;
        let (rows, pointer) = {
            let snapshot = owner.borrow();
            let coordinates = snapshot.structure.coordinates();
            (coordinates.len(), coordinates.as_ptr().cast::<f32>())
        };
        rows.checked_mul(3).ok_or_else(|| {
            pyo3::exceptions::PyOverflowError::new_err("coordinate array shape overflows usize")
        })?;
        // SAFETY: `SnapshotOwner` retains the immutable Arc-backed snapshot,
        // `[f32; 3]` is contiguous, and the resulting array is made read-only.
        let view = unsafe { ArrayView2::from_shape_ptr((rows, 3), pointer) };
        // SAFETY: the NumPy base object is exactly the owner used to obtain the
        // pointer, so Python cannot outlive or mutate the backing allocation.
        let array = unsafe { PyArray2::borrow_from_array(&view, owner.into_any()) };
        array.readwrite().make_nonwriteable();
        Ok(array)
    }
}
#[pymethods]
impl PyStructure {
    fn __repr__(&self) -> String {
        self.inner.to_string()
    }

    #[getter]
    fn coordinates<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        self.coordinates_view(py)
    }
    #[getter]
    fn atom_count(&self) -> u32 {
        self.inner.atom_count()
    }

    #[getter]
    fn residue_count(&self) -> usize {
        self.inner.residue_count()
    }

    #[getter]
    fn chain_count(&self) -> usize {
        self.inner.chain_count()
    }

    #[getter]
    fn model_count(&self) -> usize {
        self.inner.model_count()
    }

    #[getter]
    fn atoms(&self) -> PyAtoms {
        PyAtoms::all(self.clone())
    }

    #[getter]
    fn residues(&self) -> PyResidues {
        PyResidues::all(self.clone())
    }

    #[getter]
    fn chains(&self) -> PyChains {
        PyChains::all(self.clone())
    }

    #[getter]
    fn models(&self) -> PyModels {
        PyModels::new(self.clone())
    }

    #[pyo3(signature = (query, *, policy=None))]
    fn select(
        &self,
        py: Python<'_>,
        query: &Bound<'_, PyAny>,
        policy: Option<PyRef<'_, crate::policy::PyAnalysisPolicy>>,
    ) -> PyResult<PySelection> {
        let policy =
            policy.map_or_else(molframe::AnalysisPolicy::default, |policy| policy.0.clone());
        if let Ok(query) = query.extract::<PyRef<'_, PyQuery>>() {
            return crate::policy::select_compiled(py, self, &query.compiled, &policy);
        }
        let source = query.extract::<&str>()?;
        let compiled = molframe::Query::compile(source).map_err(|findings| {
            crate::query_messages::query_error(&molframe::Findings::from(findings), source)
        })?;
        crate::policy::select_compiled(py, self, &compiled, &policy)
    }
    fn edit(&self) -> PyStructureEditor {
        PyStructureEditor {
            inner: Some(self.inner.edit()),
        }
    }

    /// Covalent bonds the structure carries.
    #[getter]
    fn bond_count(&self) -> usize {
        self.inner.engine().data().bonds.len()
    }

    /// A copy with bonds inferred from covalent radii and distances.
    #[cfg(feature = "analysis")]
    #[pyo3(signature = (*, scale=molframe::DEFAULT_BOND_RADIUS_SCALE, lower_bound=molframe::DEFAULT_MINIMUM_BOND_DISTANCE, across_chains=true, context=None))]
    fn infer_bonds(
        &self,
        py: Python<'_>,
        scale: f32,
        lower_bound: f32,
        across_chains: bool,
        context: Option<&crate::execution::PyExecutionContext>,
    ) -> PyResult<Self> {
        bonds::infer(py, self, scale, lower_bound, across_chains, context)
    }

    fn _molframe_source_v2<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<Bound<'py, pyo3::types::PyCapsule>> {
        crate::native_source::capsule(py, &self.inner)
    }
    /// Per-residue states in topology order, including deposited assignments.
    #[getter]
    fn secondary_structure(&self) -> Vec<crate::secondary::PySecondaryStructure> {
        self.inner
            .secondary_structure()
            .iter()
            .copied()
            .map(Into::into)
            .collect()
    }
    fn _molframe_secondary_structure(&self) -> Vec<u8> {
        self.inner
            .engine()
            .secondary_structure()
            .iter()
            .map(|state| state.code())
            .collect()
    }
    fn _molframe_secondary_source(&self) -> Vec<u8> {
        self.inner
            .secondary_source()
            .iter()
            .map(|source| source.code())
            .collect()
    }
}
#[derive(Debug)]
#[pyclass(name = "StructureEditor", skip_from_py_object)]
pub(crate) struct PyStructureEditor {
    inner: Option<molframe::StructureEditor>,
}
#[pymethods]
impl PyStructureEditor {
    fn rename_chain(&mut self, chain: u32, label: &str) -> PyResult<()> {
        let Some(editor) = self.inner.as_mut() else {
            return Err(crate::error::internal("editor has already been finished"));
        };
        editor
            .rename_chain(molframe::ChainIndex::new(chain), label)
            .map_err(crate::error::kernel)
    }

    fn finish(&mut self) -> PyResult<PyStructure> {
        let Some(editor) = self.inner.take() else {
            return Err(crate::error::internal("editor has already been finished"));
        };
        editor
            .finish()
            .map(PyStructure::new)
            .map_err(|diagnostics| crate::error::value(format!("{diagnostics:?}")))
    }
}
#[derive(Clone, Debug)]
#[pyclass(name = "Selection", frozen, skip_from_py_object)]
pub(crate) struct PySelection {
    parent: PyStructure,
    selection: molframe::Selection,
    indices: Vec<u32>,
}
#[pymethods]
impl PySelection {
    fn __len__(&self) -> usize {
        match usize::try_from(self.selection.len()) {
            Ok(value) => value,
            Err(_) => usize::MAX,
        }
    }

    fn is_stale_for(&self, structure: &PyStructure) -> bool {
        self.selection.is_stale_for(&structure.inner)
    }

    #[getter]
    fn indices<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u32>> {
        self.indices.clone().into_pyarray(py)
    }

    #[pyo3(signature = (*, model=0))]
    fn to_coordinates<'py>(
        &self,
        py: Python<'py>,
        model: u32,
    ) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let coordinates = self
            .selection
            .to_coordinates(molframe::ModelIndex::new(model));
        let rows = coordinates.len();
        let flat = coordinates.into_iter().flatten().collect::<Vec<_>>();
        let array = flat.into_pyarray(py);
        array.reshape((rows, 3))
    }
    /// Return the residues covered by this selection.
    ///
    /// With no argument, the selection's retained parent snapshot is used.
    /// Supplying a structure keeps the explicit stale-snapshot check available
    /// for callers that apply a selection across structure versions.
    #[pyo3(signature = (structure=None))]
    fn residues(&self, structure: Option<&PyStructure>) -> PyResult<PyResidueSelection> {
        let structure = structure.map_or(&self.parent, |value| value);
        if self.selection.is_stale_for(&structure.inner) {
            return Err(crate::error::value(
                "selection is stale for the supplied structure",
            ));
        }
        Ok(PyResidueSelection::from_selection(
            structure.clone(),
            &self.selection,
        ))
    }

    fn __or__(&self, other: &Self) -> Self {
        let combined = &self.selection | &other.selection;
        Self::from_native(self.parent.clone(), combined)
    }

    fn __and__(&self, other: &Self) -> Self {
        let combined = &self.selection & &other.selection;
        Self::from_native(self.parent.clone(), combined)
    }

    fn __sub__(&self, other: &Self) -> Self {
        let combined = &self.selection - &other.selection;
        Self::from_native(self.parent.clone(), combined)
    }
}

impl PySelection {
    pub(crate) fn from_native(parent: PyStructure, selection: molframe::Selection) -> Self {
        let indices = selection.atoms().map(|atom| atom.index().get()).collect();
        Self {
            parent,
            selection,
            indices,
        }
    }
}

#[derive(Clone, Debug)]
#[pyclass(name = "Query", frozen, skip_from_py_object)]
pub(crate) struct PyQuery {
    compiled: molframe::Query,
}

impl PyQuery {
    pub(crate) fn from_native(compiled: molframe::Query) -> Self {
        Self { compiled }
    }

    pub(crate) fn native(&self) -> &molframe::Query {
        &self.compiled
    }
}

#[pymethods]
impl PyQuery {
    #[new]
    fn new(source: &str) -> PyResult<Self> {
        let compiled = molframe::Query::compile(source).map_err(|findings| {
            crate::query_messages::query_error(&molframe::Findings::from(findings), source)
        })?;
        Ok(Self { compiled })
    }

    #[pyo3(signature = (structure, *, policy=None))]
    fn select(
        &self,
        py: Python<'_>,
        structure: &PyStructure,
        policy: Option<PyRef<'_, crate::policy::PyAnalysisPolicy>>,
    ) -> PyResult<PySelection> {
        let policy =
            policy.map_or_else(molframe::AnalysisPolicy::default, |policy| policy.0.clone());
        crate::policy::select_compiled(py, structure, &self.compiled, &policy)
    }

    /// Stable identity of the normalized typed query plan.
    #[getter]
    fn fingerprint(&self) -> String {
        self.compiled.fingerprint().to_string()
    }

    /// Named queries this query refers to with `$name`, sorted.
    #[getter]
    fn references(&self) -> Vec<String> {
        self.compiled
            .references()
            .into_iter()
            .map(str::to_owned)
            .collect()
    }

    /// Canonical text accepted by every `MolFrame` query consumer.
    #[getter]
    fn source(&self) -> &str {
        self.compiled.source()
    }

    fn __and__(&self, other: &Self) -> Self {
        Self::from_native(self.compiled.clone() & other.compiled.clone())
    }

    fn __or__(&self, other: &Self) -> Self {
        Self::from_native(self.compiled.clone() | other.compiled.clone())
    }

    fn __invert__(&self) -> Self {
        Self::from_native(!self.compiled.clone())
    }
}

struct PythonBytes(PyBackedBytes);

impl AsRef<[u8]> for PythonBytes {
    fn as_ref(&self) -> &[u8] {
        self.0.as_ref()
    }
}

impl fmt::Debug for PythonBytes {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("PythonBytes")
            .field("len", &self.0.len())
            .finish()
    }
}

#[derive(Clone, Debug)]
#[pyclass(name = "Reader", frozen, skip_from_py_object)]
pub(crate) struct PyReader {
    input: molframe::InputBuffer,
    name: Option<Box<str>>,
}

#[pymethods]
impl PyReader {
    #[new]
    #[pyo3(signature = (data, *, name=None))]
    fn new(data: PyBackedBytes, name: Option<&str>) -> Self {
        Self {
            input: molframe::InputBuffer::from_owner(PythonBytes(data)),
            name: name.map(Into::into),
        }
    }

    fn read(&self, py: Python<'_>) -> PyResult<PyStructure> {
        let input = self.input.clone();
        let name = self.name.clone();
        py.detach(move || {
            molframe::read_buffer(&input, name.as_deref(), &molframe::ReadOptions::new())
        })
        .map(|(structure, _)| PyStructure::new(structure))
        .map_err(|findings| findings_error(&findings))
    }

    #[getter]
    fn byte_length(&self) -> usize {
        self.input.len()
    }
}

#[pyfunction]
#[pyo3(signature = (source, *, name=None))]
pub(crate) fn read(
    py: Python<'_>,
    source: &Bound<'_, PyAny>,
    name: Option<&str>,
) -> PyResult<PyStructure> {
    if let Ok(path) = source.extract::<PathBuf>() {
        return py
            .detach(move || molframe::read(path))
            .map(PyStructure::new)
            .map_err(|findings| findings_error(&findings));
    }
    let bytes = source.extract::<PyBackedBytes>()?;
    PyReader::new(bytes, name).read(py)
}

#[cfg(feature = "geometry")]
pub(crate) fn coordinates<'a>(array: &'a PyReadonlyArray2<'_, f32>) -> PyResult<&'a [[f32; 3]]> {
    let shape = array.shape();
    if shape.len() != 2 || shape[1] != 3 {
        return Err(crate::error::value("coordinates must have shape (n, 3)"));
    }
    let contiguous = array.as_slice().map_err(|_| {
        crate::error::value("coordinates must be C-contiguous; call numpy.ascontiguousarray")
    })?;
    Ok(contiguous.as_chunks::<3>().0)
}

#[pyfunction]
#[cfg(feature = "geometry")]
pub(crate) fn centroid(array: &Bound<'_, PyArray2<f32>>) -> PyResult<Option<[f64; 3]>> {
    let array = array.readonly();
    Ok(molframe::geometry::centroid(coordinates(&array)?))
}

#[pyfunction]
#[cfg(feature = "geometry")]
pub(crate) fn rmsd(
    mobile: &Bound<'_, PyArray2<f32>>,
    reference: &Bound<'_, PyArray2<f32>>,
) -> PyResult<f64> {
    let mobile = mobile.readonly();
    let reference = reference.readonly();
    molframe::geometry::rmsd(coordinates(&mobile)?, coordinates(&reference)?)
        .map_err(crate::error::kernel)
}

#[pyfunction]
#[cfg(feature = "geometry")]
pub(crate) fn distance_matrix<'py>(
    py: Python<'py>,
    array: &Bound<'_, PyArray2<f32>>,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let array = array.readonly();
    let matrix =
        molframe::geometry::distance_matrix(coordinates(&array)?).map_err(crate::error::kernel)?;
    let rows = matrix.rows();
    matrix.into_values().into_pyarray(py).reshape((rows, rows))
}

pub(crate) fn findings_error(findings: &molframe::Findings) -> PyErr {
    crate::error::from_findings(findings.as_slice(), &findings.to_string())
}
