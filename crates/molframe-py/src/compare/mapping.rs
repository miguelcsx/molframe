//! CCD-backed chain-sequence extraction and deterministic mapping.

use super::PyScoring;
use crate::chemistry::PyComponentDictionary;
use crate::contract::{PyAnalysis, analysis_with_value};
use crate::query::{PyAnalysisPolicy, PyNamespace};
use crate::structure::PyStructure;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyList;

#[pyclass(name = "ChainSequence", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyChainSequence {
    #[pyo3(get)]
    label: String,
    #[pyo3(get)]
    sequence: Vec<u8>,
}

#[pyclass(name = "ChainMapping", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyChainMapping {
    #[pyo3(get)]
    reference: String,
    #[pyo3(get)]
    target: String,
    #[pyo3(get)]
    identity: f64,
}

#[pyclass(name = "ChainAlternative", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyChainAlternative {
    #[pyo3(get)]
    reference: String,
    #[pyo3(get)]
    target: String,
    #[pyo3(get)]
    identity: f64,
}

#[pyclass(name = "ChainAssignment", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyChainAssignment {
    #[pyo3(get)]
    primary: Vec<PyChainMapping>,
    #[pyo3(get)]
    alternatives: Vec<PyChainAlternative>,
}

#[pyclass(name = "ResidueMatch", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyResidueMatch {
    #[pyo3(get)]
    query_position: usize,
    #[pyo3(get)]
    residue: u32,
}

#[pyfunction]
pub(crate) fn chain_sequences(
    py: Python<'_>,
    structure: &PyStructure,
    dictionary: &PyComponentDictionary,
    namespace: PyNamespace,
) -> PyResult<Vec<PyChainSequence>> {
    let structure = structure.structure().clone();
    let dictionary = dictionary.0.clone();
    py.detach(move || {
        molframe::compare::chain_sequences(&structure, dictionary.as_ref(), namespace.into())
    })
    .map(|values| values.into_iter().map(Into::into).collect())
    .map_err(value_error)
}

#[pyfunction]
pub(crate) fn map_chains(
    py: Python<'_>,
    reference: &PyStructure,
    target: &PyStructure,
    dictionary: &PyComponentDictionary,
    namespace: PyNamespace,
    scoring: PyScoring,
    minimum_identity: f64,
) -> PyResult<Vec<PyChainMapping>> {
    let reference = reference.structure().clone();
    let target = target.structure().clone();
    let dictionary = dictionary.0.clone();
    py.detach(move || {
        molframe::compare::map_chains(
            &reference,
            &target,
            dictionary.as_ref(),
            namespace.into(),
            scoring.inner(),
            minimum_identity,
        )
    })
    .map(|values| values.into_iter().map(Into::into).collect())
    .map_err(value_error)
}

#[pyfunction]
pub(crate) fn assign_chains(
    py: Python<'_>,
    reference: &PyStructure,
    target: &PyStructure,
    dictionary: &PyComponentDictionary,
    namespace: PyNamespace,
    scoring: PyScoring,
    minimum_identity: f64,
) -> PyResult<PyChainAssignment> {
    let reference = reference.structure().clone();
    let target = target.structure().clone();
    let dictionary = dictionary.0.clone();
    py.detach(move || {
        molframe::compare::assign_chains(
            &reference,
            &target,
            dictionary.as_ref(),
            namespace.into(),
            scoring.inner(),
            minimum_identity,
        )
    })
    .map(Into::into)
    .map_err(value_error)
}

#[pyfunction]
pub(crate) fn analyse_assign_chains(
    py: Python<'_>,
    reference: &PyStructure,
    target: &PyStructure,
    dictionary: &PyComponentDictionary,
    scoring: PyScoring,
    minimum_identity: f64,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let reference = reference.structure().clone();
    let target = target.structure().clone();
    let dictionary = dictionary.0.clone();
    let policy = policy.inner.clone();
    let analysis = py
        .detach(move || {
            molframe::compare::governed_assign_chains(
                &reference,
                &target,
                dictionary.as_ref(),
                scoring.inner(),
                minimum_identity,
                &policy,
            )
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, assignment_value)
}

#[pyfunction]
pub(crate) fn map_sequence_to_structure(
    py: Python<'_>,
    query: Vec<u8>,
    structure: &PyStructure,
    dictionary: &PyComponentDictionary,
    namespace: PyNamespace,
    scoring: PyScoring,
) -> PyResult<Vec<PyResidueMatch>> {
    let structure = structure.structure().clone();
    let dictionary = dictionary.0.clone();
    py.detach(move || {
        molframe::compare::map_sequence_to_structure(
            &query,
            &structure,
            dictionary.as_ref(),
            namespace.into(),
            scoring.inner(),
        )
    })
    .map(|values| values.into_iter().map(Into::into).collect())
    .map_err(value_error)
}

#[pyfunction]
pub(crate) fn analyse_map_sequence_to_structure(
    py: Python<'_>,
    query: Vec<u8>,
    structure: &PyStructure,
    dictionary: &PyComponentDictionary,
    scoring: PyScoring,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let structure = structure.structure().clone();
    let dictionary = dictionary.0.clone();
    let policy = policy.inner.clone();
    let analysis = py
        .detach(move || {
            molframe::compare::governed_map_sequence_to_structure(
                &query,
                &structure,
                dictionary.as_ref(),
                scoring.inner(),
                &policy,
            )
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, residue_match_list)
}

impl From<molframe::compare::ChainSequence> for PyChainSequence {
    fn from(value: molframe::compare::ChainSequence) -> Self {
        Self {
            label: value.label,
            sequence: value.sequence,
        }
    }
}

impl From<molframe::compare::ChainMapping> for PyChainMapping {
    fn from(value: molframe::compare::ChainMapping) -> Self {
        Self {
            reference: value.reference,
            target: value.target,
            identity: value.identity,
        }
    }
}

impl From<molframe::compare::ChainAlternative> for PyChainAlternative {
    fn from(value: molframe::compare::ChainAlternative) -> Self {
        Self {
            reference: value.reference,
            target: value.target,
            identity: value.identity,
        }
    }
}

impl From<molframe::compare::ChainAssignment> for PyChainAssignment {
    fn from(value: molframe::compare::ChainAssignment) -> Self {
        Self {
            primary: value.primary.into_iter().map(Into::into).collect(),
            alternatives: value.alternatives.into_iter().map(Into::into).collect(),
        }
    }
}

fn assignment_value(
    py: Python<'_>,
    value: molframe::compare::ChainAssignment,
) -> PyResult<Py<PyAny>> {
    Ok(Py::new(py, PyChainAssignment::from(value))?.into_any())
}

impl From<molframe::compare::ResidueMatch> for PyResidueMatch {
    fn from(value: molframe::compare::ResidueMatch) -> Self {
        Self {
            query_position: value.query_position,
            residue: value.residue.get(),
        }
    }
}

fn residue_match_list(
    py: Python<'_>,
    values: Vec<molframe::compare::ResidueMatch>,
) -> PyResult<Py<PyAny>> {
    let items = values
        .into_iter()
        .map(|value| Py::new(py, PyResidueMatch::from(value)))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(PyList::new(py, items)?.unbind().into_any())
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
