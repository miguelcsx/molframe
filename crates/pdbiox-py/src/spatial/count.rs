//! Borrowed coordinates and charged selection copies for scalar pair counts.

use super::bindings::PySpatialSearchOptions;
use crate::core::execution::PyExecutionContext;
use crate::crystallography::PyUnitCell;
use crate::geometry::borrowed_coordinates;
use numpy::{PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::{PyMemoryError, PyValueError};
use pyo3::prelude::*;

#[pyfunction]
#[pyo3(signature = (positions, cutoff, context, *, left=None, right=None, options=None, cell=None))]
fn count_pairs_within(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    cutoff: f32,
    context: &PyExecutionContext,
    left: Option<PyReadonlyArray1<'_, u32>>,
    right: Option<PyReadonlyArray1<'_, u32>>,
    options: Option<PySpatialSearchOptions>,
    cell: Option<&PyUnitCell>,
) -> PyResult<u64> {
    let positions = borrowed_coordinates(&positions)?;
    let left = left.as_ref().map(PyReadonlyArray1::as_slice).transpose()?;
    let right = right.as_ref().map(PyReadonlyArray1::as_slice).transpose()?;
    let options = options.map_or_else(
        pdbiox::SpatialSearchOptions::default,
        PySpatialSearchOptions::inner,
    );
    let periodic = cell
        .map(|cell| pdbiox::PeriodicBox::from_cell(cell.cell))
        .transpose()
        .map_err(value_error)?;
    let context = context.native();
    py.detach(move || {
        let bytes = left
            .map_or(0, <[u32]>::len)
            .checked_add(right.map_or(0, <[u32]>::len))
            .and_then(|count| count.checked_mul(size_of::<u32>()))
            .ok_or_else(|| PyMemoryError::new_err("selection size overflow"))?;
        let _reservation = context
            .try_reserve(bytes)
            .map_err(|error| PyMemoryError::new_err(error.to_string()))?;
        let atoms = u32::try_from(positions.len()).map_err(value_error)?;
        let selection = |indices: Option<&[u32]>| {
            indices.map_or(pdbiox::AtomSelection::All(atoms), |indices| {
                pdbiox::AtomSelection::Sparse(indices.to_vec())
            })
        };
        let left = selection(left);
        let right = selection(right);
        pdbiox::spatial::count_pairs_within(&pdbiox::spatial::PairQuery {
            positions,
            left: &left,
            right: &right,
            cutoff,
            options,
            periodic: periodic.as_ref(),
            context: &context,
        })
        .map_err(|error| match error {
            pdbiox::SpatialError::Memory(_) => PyMemoryError::new_err(error.to_string()),
            _ => value_error(error),
        })
    })
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(count_pairs_within, module)?)
}
