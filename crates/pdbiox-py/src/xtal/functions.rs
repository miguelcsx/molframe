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
    pdbiox::xtal::lower_assemblies(&document.inner)
        .map(PyAssemblySet)
        .map_err(|findings| read_error(py, &findings))
}

#[pyfunction]
pub(crate) fn lower_symmetry(py: Python<'_>, document: &PyCifDocument) -> PyResult<PySymmetrySet> {
    pdbiox::xtal::lower_symmetry(&document.inner)
        .map(PySymmetrySet)
        .map_err(|findings| read_error(py, &findings))
}

#[pyfunction]
#[pyo3(signature = (structure, symmetry, model=0, cutoff=5.0, *, limit=None, backend=PySpatialBackend::Auto))]
pub(crate) fn crystal_neighbors(
    structure: &PyStructure,
    symmetry: &PySymmetrySet,
    model: usize,
    cutoff: f64,
    limit: Option<usize>,
    backend: PySpatialBackend,
) -> PyResult<Vec<PyCrystalNeighbor>> {
    let model = model_index(model)?;
    let result = match limit {
        Some(limit) => pdbiox::xtal::crystal_neighbors_with_limit(
            structure.structure(),
            &symmetry.0,
            model,
            cutoff,
            limit,
        ),
        None => pdbiox::xtal::crystal_neighbors_with_backend(
            structure.structure(),
            &symmetry.0,
            model,
            cutoff,
            backend.into(),
            pdbiox::xtal::DEFAULT_CRYSTAL_IMAGE_LIMIT,
        ),
    };
    result
        .map(|values| values.into_iter().map(Into::into).collect())
        .map_err(value_error)
}

#[pyfunction]
#[pyo3(signature = (structure, symmetry, model=0, cutoff=5.0, *, limit))]
pub(crate) fn crystal_neighbors_with_limit(
    structure: &PyStructure,
    symmetry: &PySymmetrySet,
    model: usize,
    cutoff: f64,
    limit: usize,
) -> PyResult<Vec<PyCrystalNeighbor>> {
    crystal_neighbors(
        structure,
        symmetry,
        model,
        cutoff,
        Some(limit),
        PySpatialBackend::Auto,
    )
}

#[pyfunction]
#[pyo3(signature = (structure, symmetry, model=0, cutoff=5.0, backend=PySpatialBackend::Auto))]
pub(crate) fn crystal_neighbors_with_backend(
    structure: &PyStructure,
    symmetry: &PySymmetrySet,
    model: usize,
    cutoff: f64,
    backend: PySpatialBackend,
) -> PyResult<Vec<PyCrystalNeighbor>> {
    crystal_neighbors(structure, symmetry, model, cutoff, None, backend)
}

fn model_index(model: usize) -> PyResult<pdbiox::ModelIndex> {
    u32::try_from(model)
        .map(pdbiox::ModelIndex::new)
        .map_err(|_| value_error("model index exceeds u32"))
}

fn value_error(error: impl std::fmt::Display) -> pyo3::PyErr {
    pyo3::exceptions::PyValueError::new_err(error.to_string())
}
