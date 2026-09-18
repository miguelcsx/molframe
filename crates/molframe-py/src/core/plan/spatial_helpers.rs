//! Input, option and serialization helpers for spatial operation nodes.

use crate::crystallography::PyUnitCell;
use crate::geometry::borrowed_coordinates;
use crate::spatial::{PySpatialSearchOptions, spatial_arrays::validate_sorted_indices};
use numpy::{PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::{PyKeyError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};
use std::collections::HashMap;

pub(super) fn position_count(py: Python<'_>, positions: &Py<PyAny>) -> PyResult<usize> {
    let array = positions.bind(py).extract::<PyReadonlyArray2<'_, f32>>()?;
    let coordinates = borrowed_coordinates(&array)?;
    Ok(coordinates.len())
}

pub(super) fn validate_optional_indices(
    py: Python<'_>,
    value: Option<&Py<PyAny>>,
    count: usize,
) -> PyResult<()> {
    if let Some(value) = value {
        validate_indices(py, value, count)?;
    }
    Ok(())
}

pub(super) fn validate_indices(py: Python<'_>, value: &Py<PyAny>, count: usize) -> PyResult<()> {
    let array = value.bind(py).extract::<PyReadonlyArray1<'_, u32>>()?;
    let indices = contiguous_indices(&array)?;
    validate_sorted_indices(indices, count).map_err(PyValueError::new_err)
}

/// Validates the search radius before the operation enters a reusable plan.
pub(super) fn validate_cutoff(cutoff: f32) -> PyResult<()> {
    if cutoff.is_finite() && cutoff >= 0.0 {
        Ok(())
    } else {
        Err(PyValueError::new_err(
            "cutoff must be finite and non-negative",
        ))
    }
}

/// Materializes the one owned selection required by the native reusable plan.
///
/// The input array itself remains borrowed through execution. `SpatialRequest`
/// owns an adaptive `AtomSelection`, so this bounded O(k) copy is the explicit
/// compilation boundary that lets the facade cache and reuse its index safely.
pub(super) fn owned_selection(
    py: Python<'_>,
    value: Option<&Py<PyAny>>,
    count: usize,
) -> PyResult<molframe::AtomSelection> {
    let Some(value) = value else {
        return u32::try_from(count)
            .map(molframe::AtomSelection::All)
            .map_err(|_| PyValueError::new_err("coordinate count exceeds u32 capacity"));
    };
    let array = value.bind(py).extract::<PyReadonlyArray1<'_, u32>>()?;
    let indices = contiguous_indices(&array)?;
    py.detach(move || {
        validate_sorted_indices(indices, count).map_err(PyValueError::new_err)?;
        Ok(molframe::AtomSelection::from_sorted(indices.to_vec()))
    })
}

fn contiguous_indices<'py>(array: &'py PyReadonlyArray1<'py, u32>) -> PyResult<&'py [u32]> {
    array.as_slice().map_err(|_| {
        PyValueError::new_err(
            "atom indices must be C-contiguous uint32 arrays; pass a contiguous array explicitly",
        )
    })
}

pub(super) fn spatial_options(
    py: Python<'_>,
    options: Option<&Py<PyAny>>,
) -> PyResult<PySpatialSearchOptions> {
    let Some(options) = options else {
        return Ok(PySpatialSearchOptions::balanced_value());
    };
    Ok(options
        .bind(py)
        .extract::<PyRef<'_, PySpatialSearchOptions>>()
        .map(|value| *value)?)
}

pub(super) fn validate_cell(py: Python<'_>, cell: Option<&Py<PyAny>>) -> PyResult<()> {
    if let Some(cell) = cell {
        cell.bind(py).extract::<PyRef<'_, PyUnitCell>>()?;
    }
    Ok(())
}

pub(super) fn periodic_box(
    py: Python<'_>,
    cell: Option<&Py<PyAny>>,
) -> PyResult<Option<molframe::PeriodicBox>> {
    let Some(cell) = cell else {
        return Ok(None);
    };
    let cell = cell.bind(py).extract::<PyRef<'_, PyUnitCell>>()?;
    molframe::PeriodicBox::from_cell(cell.cell)
        .map(Some)
        .map_err(|error| PyValueError::new_err(error.to_string()))
}

pub(super) fn retain_positions<'py>(
    py: Python<'py>,
    value: &Py<PyAny>,
    arrays: &mut Vec<PyReadonlyArray2<'py, f32>>,
    coordinate_slots: &mut HashMap<usize, usize>,
) -> PyResult<usize> {
    let key = value.bind(py).as_ptr() as usize;
    if let Some(slot) = coordinate_slots.get(&key).copied() {
        return Ok(slot);
    }
    let array = value.bind(py).extract::<PyReadonlyArray2<'py, f32>>()?;
    borrowed_coordinates(&array)?;
    let slot = arrays.len();
    arrays.push(array);
    coordinate_slots.insert(key, slot);
    Ok(slot)
}

pub(super) fn optional_object(
    config: &Bound<'_, PyDict>,
    name: &str,
) -> PyResult<Option<Py<PyAny>>> {
    match config.get_item(name)? {
        Some(value) if !value.is_none() => Ok(Some(value.unbind())),
        _ => Ok(None),
    }
}

pub(super) fn set_optional_object(
    result: &Bound<'_, PyDict>,
    name: &str,
    py: Python<'_>,
    value: Option<&Py<PyAny>>,
) -> PyResult<()> {
    match value {
        Some(value) => result.set_item(name, value.bind(py)),
        None => result.set_item(name, py.None()),
    }
}

pub(super) fn required<'py>(
    config: &'py Bound<'py, PyDict>,
    name: &str,
) -> PyResult<Bound<'py, PyAny>> {
    config
        .get_item(name)?
        .ok_or_else(|| PyKeyError::new_err(format!("missing operation field: {name}")))
}

pub(super) fn incompatible_result() -> PyErr {
    PyValueError::new_err("native spatial plan returned an incompatible result")
}
