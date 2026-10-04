//! The chain of one structure that corresponds to each chain of another, and the residues of
//! a structure that a query sequence lands on, both by sequence alignment.

use crate::bindings::PyStructure;
use crate::sequence::PyScoring;
use molframe::Namespace;
use pyo3::prelude::*;
use std::path::PathBuf;

fn mapping_inputs(
    components: &std::path::Path,
    version: &str,
) -> PyResult<molframe::chemistry::CifProvider> {
    molframe::read_component_dictionary(components, molframe::DictionaryVersion::new(version))
        .map(|(provider, _)| provider)
        .map_err(|findings| crate::bindings::findings_error(&findings))
}

/// The chain of one structure that each chain of another corresponds to, by sequence.
///
/// Returns `(reference, target, identity)` for the assignment of the greatest total identity
/// and, under `alternatives`, every other pair that meets `min_identity`.
#[pyfunction]
#[pyo3(signature = (
    reference,
    target,
    *,
    components,
    components_version,
    scoring,
    min_identity,
    chain_names="label",
))]
#[allow(clippy::too_many_arguments, clippy::needless_pass_by_value)]
fn assign_chains(
    py: Python<'_>,
    reference: &PyStructure,
    target: &PyStructure,
    components: PathBuf,
    components_version: &str,
    scoring: &PyScoring,
    min_identity: f64,
    chain_names: &str,
) -> PyResult<PyChainAssignment> {
    let provider = mapping_inputs(&components, components_version)?;
    let namespace: Namespace = chain_names.parse().map_err(crate::error::kernel)?;
    let (reference, target) = (reference.inner.clone(), target.inner.clone());
    let scoring = scoring.0;
    let found = py
        .detach(|| {
            molframe::compare::assign_chains(
                reference.engine(),
                target.engine(),
                &provider,
                namespace,
                scoring,
                min_identity,
            )
        })
        .map_err(crate::error::kernel)?;
    Ok(PyChainAssignment { inner: found })
}

/// A one-to-one chain assignment and the alternatives to it.
#[derive(Clone, Debug)]
#[pyclass(
    name = "ChainAssignment",
    frozen,
    skip_from_py_object,
    module = "molframe.compare"
)]
pub(crate) struct PyChainAssignment {
    inner: molframe::compare::ChainAssignment,
}

#[pymethods]
impl PyChainAssignment {
    /// `(reference, target, identity)` of the assignment of the greatest total identity.
    #[getter]
    fn primary(&self) -> Vec<(String, String, f64)> {
        self.inner
            .primary
            .iter()
            .map(|pair| (pair.reference.clone(), pair.target.clone(), pair.identity))
            .collect()
    }

    /// `(reference, target, identity)` of the pairs left out that still meet the threshold.
    #[getter]
    fn alternatives(&self) -> Vec<(String, String, f64)> {
        self.inner
            .alternatives
            .iter()
            .map(|pair| (pair.reference.clone(), pair.target.clone(), pair.identity))
            .collect()
    }
}

/// Where each residue of a query sequence lies in a structure, by alignment.
///
/// Returns `(query_position, residue_index)` pairs.
#[pyfunction]
#[pyo3(signature = (
    query,
    structure,
    *,
    components,
    components_version,
    scoring,
    chain_names="label",
))]
#[allow(clippy::needless_pass_by_value)]
fn map_sequence_to_structure(
    py: Python<'_>,
    query: &str,
    structure: &PyStructure,
    components: PathBuf,
    components_version: &str,
    scoring: &PyScoring,
    chain_names: &str,
) -> PyResult<Vec<(usize, u32)>> {
    let provider = mapping_inputs(&components, components_version)?;
    let namespace: Namespace = chain_names.parse().map_err(crate::error::kernel)?;
    let structure = structure.inner.clone();
    let scoring = scoring.0;
    let query = query.as_bytes().to_vec();
    let found = py
        .detach(|| {
            molframe::compare::map_sequence_to_structure(
                &query,
                structure.engine(),
                &provider,
                namespace,
                scoring,
            )
        })
        .map_err(crate::error::kernel)?;
    Ok(found
        .into_iter()
        .map(|found| (found.query_position, found.residue.get()))
        .collect())
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyChainAssignment>()?;
    module.add_function(wrap_pyfunction!(assign_chains, module)?)?;
    module.add_function(wrap_pyfunction!(map_sequence_to_structure, module)?)
}
