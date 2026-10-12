use super::PyStructure;
use crate::hierarchy::PyResidueSelection;
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
