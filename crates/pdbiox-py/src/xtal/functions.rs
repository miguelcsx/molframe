//! Document lowering and crystal-neighbour entry points.

use super::assembly::PyAssemblySet;
use super::ncs::PyCrystalNeighbor;
use super::types::PySymmetrySet;
use crate::cif_document::PyCifDocument;
use crate::errors::read_error;
use crate::graph::PySpatialBackend;
use crate::structure::PyStructure;
use pyo3::prelude::*;

#[pyfunction]
pub(crate) fn lower_assemblies(
    py: Python<'_>,
    document: &PyCifDocument,
) -> PyResult<PyAssemblySet> {
    // Only the diagnostic conversion needs the interpreter.
    py.detach(|| pdbiox::xtal::lower_assemblies(&document.inner))
        .map(PyAssemblySet)
        .map_err(|findings| read_error(py, &findings))
}

#[pyfunction]
pub(crate) fn lower_symmetry(py: Python<'_>, document: &PyCifDocument) -> PyResult<PySymmetrySet> {
    // Only the diagnostic conversion needs the interpreter, so the lowering
    // itself runs with the lock released.
    py.detach(|| pdbiox::xtal::lower_symmetry(&document.inner))
        .map(PySymmetrySet)
        .map_err(|findings| read_error(py, &findings))
}

#[pyfunction]
#[pyo3(signature = (structure, symmetry, model=0, cutoff=5.0, *, limit=None, backend=PySpatialBackend::Auto))]
pub(crate) fn collect_crystal_neighbors(
    py: Python<'_>,
    structure: &PyStructure,
    symmetry: &PySymmetrySet,
    model: usize,
    cutoff: f64,
    limit: Option<usize>,
    backend: PySpatialBackend,
) -> PyResult<Vec<PyCrystalNeighbor>> {
    py.detach(move || neighbors_of(structure, symmetry, model, cutoff, limit, backend))
}

/// Native crystal-neighbour search kept outside the interpreter attachment.
fn neighbors_of(
    structure: &PyStructure,
    symmetry: &PySymmetrySet,
    model: usize,
    cutoff: f64,
    limit: Option<usize>,
    backend: PySpatialBackend,
) -> PyResult<Vec<PyCrystalNeighbor>> {
    let model = model_index(model)?;
    let candidate_limit = match limit {
        Some(limit) => limit,
        None => pdbiox::xtal::DEFAULT_CRYSTAL_IMAGE_LIMIT,
    };
    let options = pdbiox::xtal::CrystalNeighborOptions {
        backend: backend.into(),
        candidate_limit,
        ..pdbiox::xtal::CrystalNeighborOptions::default()
    };
    let result = pdbiox::xtal::collect_crystal_neighbors(
        structure.structure(),
        &symmetry.0,
        model,
        cutoff,
        options,
        &crate::core::execution::default_context(),
    );
    result
        .map(|values| values.into_iter().map(Into::into).collect())
        .map_err(value_error)
}

fn model_index(model: usize) -> PyResult<pdbiox::ModelIndex> {
    u32::try_from(model)
        .map(pdbiox::ModelIndex::new)
        .map_err(|_| value_error("model index exceeds u32"))
}

fn value_error(error: impl std::fmt::Display) -> pyo3::PyErr {
    pyo3::exceptions::PyValueError::new_err(error.to_string())
}
