//! Remaining policy-bound comparison workflows delegated to Rust.

use super::{PyCeAlignment, PyCeOptions, PyChainMapping, PyScoring};
use crate::chemistry::PyComponentDictionary;
use crate::contract::{PyAnalysis, analysis_with_value};
use crate::geometry::borrowed_coordinates;
use crate::query::PyAnalysisPolicy;
use crate::structure::PyStructure;
use numpy::PyReadonlyArray2;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyFloat, PyList};

#[pyfunction]
pub(crate) fn analyse_gdt(
    py: Python<'_>,
    model: PyReadonlyArray2<'_, f32>,
    reference: PyReadonlyArray2<'_, f32>,
    cutoffs: Vec<f64>,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let model = borrowed_coordinates(&model)?;
    let reference = borrowed_coordinates(&reference)?;
    let cutoffs = cutoffs.into_boxed_slice();
    let policy = policy.inner.clone();
    let analysis = py
        .detach(|| {
            molframe::compare::governed_gdt_with_cutoffs(model, reference, &cutoffs, &policy)
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, |py, value| {
        Ok(PyFloat::new(py, value).unbind().into_any())
    })
}

#[pyfunction]
pub(crate) fn analyse_weighted_rmsd(
    py: Python<'_>,
    model: PyReadonlyArray2<'_, f32>,
    reference: PyReadonlyArray2<'_, f32>,
    weights: Vec<f64>,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let model = borrowed_coordinates(&model)?;
    let reference = borrowed_coordinates(&reference)?;
    let policy = policy.inner.clone();
    let analysis = py
        .detach(|| molframe::compare::governed_weighted_rmsd(model, reference, &weights, &policy))
        .map_err(value_error)?;
    analysis_with_value(py, analysis, |py, value| {
        Ok(PyFloat::new(py, value).unbind().into_any())
    })
}

#[pyfunction]
pub(crate) fn analyse_ce_align(
    py: Python<'_>,
    reference: PyReadonlyArray2<'_, f32>,
    mobile: PyReadonlyArray2<'_, f32>,
    options: PyCeOptions,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let reference = borrowed_coordinates(&reference)?;
    let mobile = borrowed_coordinates(&mobile)?;
    let policy = policy.inner.clone();
    let analysis = py
        .detach(|| {
            molframe::compare::governed_ce_align(reference, mobile, options.inner(), &policy)
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, |py, value| {
        Py::new(py, PyCeAlignment::from(value)).map(Py::into_any)
    })
}

#[pyfunction]
pub(crate) fn analyse_ce_alignments(
    py: Python<'_>,
    reference: PyReadonlyArray2<'_, f32>,
    mobile: PyReadonlyArray2<'_, f32>,
    options: PyCeOptions,
    policy: &PyAnalysisPolicy,
) -> PyResult<PyAnalysis> {
    let reference = borrowed_coordinates(&reference)?;
    let mobile = borrowed_coordinates(&mobile)?;
    let policy = policy.inner.clone();
    let analysis = py
        .detach(|| {
            molframe::compare::governed_ce_alignments(reference, mobile, options.inner(), &policy)
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, ce_alignment_list)
}

#[pyfunction]
pub(crate) fn analyse_map_chains(
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
            molframe::compare::governed_map_chains(
                &reference,
                &target,
                dictionary.as_ref(),
                scoring.inner(),
                minimum_identity,
                &policy,
            )
        })
        .map_err(value_error)?;
    analysis_with_value(py, analysis, chain_mapping_list)
}

fn ce_alignment_list(
    py: Python<'_>,
    values: Vec<molframe::compare::CeAlignment>,
) -> PyResult<Py<PyAny>> {
    let items = values
        .into_iter()
        .map(|value| Py::new(py, PyCeAlignment::from(value)))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(PyList::new(py, items)?.unbind().into_any())
}

fn chain_mapping_list(
    py: Python<'_>,
    values: Vec<molframe::compare::ChainMapping>,
) -> PyResult<Py<PyAny>> {
    let items = values
        .into_iter()
        .map(|value| Py::new(py, PyChainMapping::from(value)))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(PyList::new(py, items)?.unbind().into_any())
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
