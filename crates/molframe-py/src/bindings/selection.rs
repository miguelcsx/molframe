//! Python wrappers for atom selections and compiled queries.
//!
//! Complexity: O(k) time and space for k selected atoms when exporting indices or
//! coordinates; set operations are linear in the operands; compiling a query is
//! linear in its source text.

use super::PyStructure;
use crate::bindings::array::Array1;
use crate::hierarchy::PyResidueSelection;
use molframe::engine::core::AtomSelection;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods};
use pyo3::prelude::*;

#[derive(Clone, Debug)]
#[pyclass(name = "Selection", frozen, skip_from_py_object)]
pub(crate) struct PySelection {
    pub(super) parent: PyStructure,
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
        let coordinates = py.detach(|| {
            self.selection
                .to_coordinates(molframe::ModelIndex::new(model))
        });
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

    /// The atoms at these positions, in any order; a repeated position counts once.
    #[staticmethod]
    fn from_indices(structure: &PyStructure, indices: Array1<'_, '_, u32>) -> PyResult<Self> {
        let count = structure.inner.engine().atom_count();
        let readonly = indices.readonly();
        let mut atoms = readonly
            .as_slice()
            .map_err(|_| crate::error::value("indices must be a contiguous uint32 array"))?
            .to_vec();
        atoms.sort_unstable();
        atoms.dedup();
        if atoms.last().is_some_and(|&last| last >= count) {
            return Err(crate::error::value(
                "an index is outside the structure's atoms",
            ));
        }
        Ok(Self::of_atoms(structure, AtomSelection::from_sorted(atoms)))
    }

    /// The atoms whose entry in `mask` is true; the mask has one entry per atom.
    #[staticmethod]
    fn from_mask(structure: &PyStructure, mask: Array1<'_, '_, bool>) -> PyResult<Self> {
        let readonly = mask.readonly();
        let flags = readonly
            .as_slice()
            .map_err(|_| crate::error::value("mask must be a contiguous boolean array"))?;
        if usize::try_from(structure.inner.engine().atom_count()) != Ok(flags.len()) {
            return Err(crate::error::value(
                "mask length must equal the structure's atom count",
            ));
        }
        let mut atoms = Vec::new();
        for (index, &selected) in flags.iter().enumerate() {
            if selected {
                let atom = u32::try_from(index)
                    .map_err(|_| crate::error::value("mask is longer than u32 atoms"))?;
                atoms.push(atom);
            }
        }
        Ok(Self::of_atoms(structure, AtomSelection::from_sorted(atoms)))
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

    /// A selection over `atoms` of `structure`, a zero-copy view of its atoms.
    fn of_atoms(structure: &PyStructure, atoms: AtomSelection) -> Self {
        let view = structure.inner.engine().view_of(atoms);
        Self::from_native(structure.clone(), molframe::Selection::from(view))
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

    #[pyo3(signature = (structure, *, policy=None, context=None))]
    fn select(
        &self,
        py: Python<'_>,
        structure: &PyStructure,
        policy: Option<PyRef<'_, crate::policy::PyAnalysisPolicy>>,
        context: Option<&crate::execution::PyExecutionContext>,
    ) -> PyResult<PySelection> {
        let policy =
            policy.map_or_else(molframe::AnalysisPolicy::default, |policy| policy.0.clone());
        crate::policy::select_compiled_in(py, structure, &self.compiled, &policy, context)
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
