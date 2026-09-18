//! Chemistry-aware motif mapping bindings.

use super::errors::mapping_error;
use super::specification::{PyAtomSite, PyMotif};
use crate::chemistry::PyComponentDictionary;
use crate::query::PyAnalysisPolicy;
use crate::structure::PyStructure;
use pyo3::prelude::*;
use std::collections::BTreeMap;

#[pyclass(name = "MappedMotif", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyMappedMotif(pub(crate) molframe::fx::MappedMotif);

#[pymethods]
impl PyMappedMotif {
    #[new]
    fn new(components: BTreeMap<String, u32>, atoms: Vec<(PyAtomSite, Vec<u32>)>) -> Self {
        let components = components
            .into_iter()
            .map(|(name, residue)| (name.into_boxed_str(), molframe::ResidueIndex::new(residue)))
            .collect();
        let atoms = atoms
            .into_iter()
            .map(|(site, values)| {
                (
                    site.0,
                    values.into_iter().map(molframe::AtomIndex::new).collect(),
                )
            })
            .collect();
        Self(molframe::fx::MappedMotif { components, atoms })
    }

    #[getter]
    fn components(&self) -> BTreeMap<String, u32> {
        self.0
            .components
            .iter()
            .map(|(name, residue)| (name.to_string(), residue.get()))
            .collect()
    }

    #[getter]
    fn atoms(&self) -> Vec<(PyAtomSite, Vec<u32>)> {
        self.0
            .atoms
            .iter()
            .map(|(site, values)| {
                (
                    PyAtomSite(site.clone()),
                    values.iter().map(|atom| atom.get()).collect(),
                )
            })
            .collect()
    }

    fn atom_indices(&self, site: &PyAtomSite) -> Option<Vec<u32>> {
        self.0
            .atoms
            .get(&site.0)
            .map(|values| values.iter().map(|atom| atom.get()).collect())
    }
}

#[pyclass(name = "MappingSet", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyMappingSet(pub(crate) molframe::fx::MappingSet);

#[pymethods]
impl PyMappingSet {
    #[new]
    fn new(mappings: Vec<PyMappedMotif>, ambiguous: bool) -> Self {
        Self(molframe::fx::MappingSet {
            mappings: mappings.into_iter().map(|mapping| mapping.0).collect(),
            ambiguous,
        })
    }

    #[getter]
    fn mappings(&self) -> Vec<PyMappedMotif> {
        self.0.mappings.iter().cloned().map(PyMappedMotif).collect()
    }

    #[getter]
    const fn ambiguous(&self) -> bool {
        self.0.ambiguous
    }
}

#[pyfunction]
#[pyo3(signature = (structure, motif, limit, *, dictionary=None, policy=None))]
pub(crate) fn map_motif(
    py: Python<'_>,
    structure: &PyStructure,
    motif: &PyMotif,
    limit: usize,
    dictionary: Option<&PyComponentDictionary>,
    policy: Option<&PyAnalysisPolicy>,
) -> PyResult<PyMappingSet> {
    let structure = structure.structure().clone();
    let motif = motif.0.clone();
    let dictionary = dictionary.map(|value| value.0.clone());
    let policy = policy.map_or_else(molframe::AnalysisPolicy::default, |value| {
        value.inner.clone()
    });
    py.detach(move || {
        let provider = dictionary
            .as_ref()
            .map(|value| value.as_ref() as &dyn molframe::ComponentProvider);
        molframe::fx::map_motif(&structure, &motif, provider, &policy, limit)
            .map(PyMappingSet)
            .map_err(mapping_error)
    })
}
