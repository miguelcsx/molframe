//! Central facade handles and stateless namespace functions.

use crate::hierarchy::{PyAtoms, PyChains, PyModels, PyResidueSelection, PyResidues};
use numpy::ndarray::ArrayView2;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods, ToPyArray};
use pyo3::prelude::*;
use pyo3::types::PyAny;

#[cfg(feature = "analysis")]
mod analysis;
#[cfg(feature = "geometry")]
mod geometry;
#[cfg(feature = "geometry")]
pub(crate) mod measures;
#[cfg(feature = "geometry")]
pub(crate) use geometry::{centroid, coordinates, distance_matrix, rmsd};
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
    /// Concatenates snapshots in input row order, preserving bonds and annotations.
    #[staticmethod]
    fn merge(py: Python<'_>, structures: Vec<PyRef<'_, Self>>) -> PyResult<Self> {
        let sources: Vec<_> = structures
            .into_iter()
            .map(|value| value.inner.clone())
            .collect();
        py.detach(|| molframe::Structure::merge(&sources))
            .map(Self::new)
            .map_err(|findings| findings_error(&findings))
    }

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
    /// Resolves alternate conformations without changing stored atom rows.
    #[pyo3(signature = (*, policy=None))]
    fn resolve_altlocs(
        &self,
        py: Python<'_>,
        policy: Option<PyRef<'_, crate::policy::PyAnalysisPolicy>>,
    ) -> PyResult<crate::analysis_result::PyAnalysis> {
        let policy = match policy {
            Some(policy) => policy.0.clone(),
            None => molframe::AnalysisPolicy::default(),
        };
        let analysis = self.inner.resolve_altlocs(&policy);
        let value = match analysis.value() {
            Some(selection) => {
                let view = self.inner.engine().view_of(selection.clone());
                let selected =
                    PySelection::from_native(self.clone(), molframe::Selection::from(view));
                Some(Py::new(py, selected)?.into_any())
            }
            None => None,
        };
        Ok(crate::analysis_result::PyAnalysis::new(&analysis, value))
    }

    /// Starts an edit: topology changes are staged and published together by
    /// `finish()`, coordinate changes go through `coordinates()`.
    fn edit(&self) -> crate::editing::PyStructureEditor {
        crate::editing::PyStructureEditor::new(self.inner.edit())
    }

    /// The covalent bonds the structure carries, as an Arrow-readable table.
    #[getter]
    fn bonds(&self) -> crate::interop::PyBondTable {
        crate::interop::PyBondTable::new(self.clone())
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
    /// Entry-level facts the file states: identifier, title, method, resolution.
    #[getter]
    fn metadata(&self) -> crate::structure_data::PyEntryMetadata {
        crate::structure_data::PyEntryMetadata {
            parent: self.clone(),
        }
    }

    /// The crystallographic cell, or `None` when the file defines none.
    ///
    /// The placeholder unit cube some files write in place of a cell is not
    /// reported as one.
    #[getter]
    fn cell(&self) -> PyResult<Option<crate::crystal::PyUnitCell>> {
        match self.inner.cell() {
            Some(cell) if !cell.is_placeholder() => {
                crate::crystal::PyUnitCell::from_cell(&cell).map(Some)
            }
            _ => Ok(None),
        }
    }

    /// Anisotropic displacement parameters, one row per atom that has them.
    ///
    /// Columns are `atom` and the `Å²` tensor `u11`, `u22`, `u33`, `u12`,
    /// `u13`, `u23`. `None` means the file did not say whether any exist; an
    /// empty table means it said there are none.
    #[getter]
    fn anisotropy(&self, py: Python<'_>) -> Option<crate::table::PyTable> {
        crate::structure_data::anisotropy(py, &self.inner)
    }

    /// The distinct species the structure contains.
    #[getter]
    fn entities(&self) -> crate::entities::PyEntities {
        crate::entities::PyEntities::new(self.clone())
    }

    /// Custom per-atom columns, such as partial charges or confidence scores.
    #[getter]
    fn annotations(&self) -> crate::structure_data::PyAnnotations {
        crate::structure_data::PyAnnotations {
            parent: self.clone(),
        }
    }

    /// Atomic number of each atom (`0` for an unknown element).
    #[getter]
    fn elements<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u8>> {
        crate::structure_data::elements(&self.inner).to_pyarray(py)
    }

    /// Occupancy of each atom; `NaN` where the file records none.
    #[getter]
    fn occupancies<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f32>> {
        crate::structure_data::occupancies(&self.inner).to_pyarray(py)
    }

    /// Temperature factor of each atom; `NaN` where the file records none.
    #[getter]
    fn b_factors<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f32>> {
        crate::structure_data::b_factors(&self.inner).to_pyarray(py)
    }

    /// Which assigner produced each residue's secondary state.
    #[getter]
    fn secondary_source(&self) -> Vec<crate::secondary::PySecondarySource> {
        self.inner
            .secondary_source()
            .iter()
            .copied()
            .map(Into::into)
            .collect()
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
    pub(crate) fn atom_indices(&self) -> &[u32] {
        &self.indices
    }

    pub(crate) const fn native(&self) -> &molframe::Selection {
        &self.selection
    }

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

pub(crate) fn findings_error(findings: &molframe::Findings) -> PyErr {
    crate::error::from_findings(findings.as_slice(), &findings.to_string())
}
