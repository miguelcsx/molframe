//! Shared binding-side projection and coordinate conversion.

use crate::query::PyNamespace;
use crate::structure::PyStructure;
use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray2};
use pdbiox::ModelIndex;
use pdbiox::adapters::{TopologyExport, TopologyExportError};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

pub(super) fn project(
    structure: &PyStructure,
    model: usize,
    namespace: PyNamespace,
) -> PyResult<TopologyExport> {
    let model = u32::try_from(model)
        .map(ModelIndex::new)
        .map_err(|_| PyValueError::new_err("model index exceeds the native index range"))?;
    TopologyExport::from_model(structure.structure(), model, namespace.into())
        .map_err(|error| projection_error(&error))
}

pub(super) fn coordinate_array<'py>(
    py: Python<'py>,
    export: &TopologyExport,
) -> PyResult<Bound<'py, PyArray2<f32>>> {
    const CARTESIAN_DIMENSIONS: usize = 3;
    let coordinates: Vec<_> = export.atoms.iter().flat_map(|atom| atom.position).collect();
    let array = Array2::from_shape_vec((export.atoms.len(), CARTESIAN_DIMENSIONS), coordinates)
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
    Ok(array.into_pyarray(py))
}

fn projection_error(error: &TopologyExportError) -> PyErr {
    PyValueError::new_err(error.to_string())
}
