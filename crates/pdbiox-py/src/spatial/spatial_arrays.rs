//! NumPy-first spatial queries with borrowed inputs and columnar results.
//!
//! Coordinates are borrowed from C-contiguous `float32` `(n, 3)` arrays and
//! selections from sorted, unique C-contiguous `uint32` arrays. Native search
//! and result-column construction run without the GIL. Pair results stay in
//! Rust-owned columns and are projected as immutable zero-copy `NumPy` views.

use super::{PyKdPeriodicOptions, PySpatialSearchOptions, normalize_selection};
use crate::crystallography::PyUnitCell;
use crate::geometry::borrowed_coordinates;
use crate::graph::PySpatialBackend;
use numpy::ndarray::{Array1, ArrayView1, ArrayView2};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::{PyMemoryError, PyValueError};
use pyo3::prelude::*;

/// Immutable, columnar neighbour rows.
///
/// Every `indices[i]` row is `[source, target]` and `distance[i]` is its
/// Euclidean separation. The same layout represents both all-pair searches and
/// nearest-neighbour searches; in the latter, every source is the query atom.
#[pyclass(name = "NeighborTable", frozen, skip_from_py_object)]
#[derive(Debug)]
pub(crate) struct PyNeighborTable {
    indices: Vec<u32>,
    distance: Vec<f32>,
}

#[pymethods]
impl PyNeighborTable {
    /// A read-only zero-copy `uint32[m, 2]` view of source and target indices.
    #[getter]
    fn indices<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyArray2<u32>>> {
        let (rows, pointer) = {
            let table = slf.borrow();
            (table.rows()?, table.indices.as_ptr())
        };
        Ok(readonly_indices(slf, rows, pointer))
    }

    /// A read-only zero-copy `float32[m]` view aligned with `indices`.
    #[getter]
    fn distance<'py>(slf: &Bound<'py, Self>) -> Bound<'py, PyArray1<f32>> {
        let pointer = {
            let table = slf.borrow();
            table.distance.as_ptr()
        };
        readonly_distances(slf, pointer)
    }

    fn __len__(&self) -> usize {
        self.distance.len()
    }

    fn __repr__(&self) -> String {
        format!("NeighborTable(rows={})", self.distance.len())
    }
}

impl PyNeighborTable {
    fn rows(&self) -> PyResult<usize> {
        let Some(expected_indices) = self.distance.len().checked_mul(2) else {
            return Err(PyValueError::new_err(
                "neighbor-table row count exceeds addressable memory",
            ));
        };
        if self.indices.len() != expected_indices {
            return Err(PyValueError::new_err(
                "neighbor-table columns have inconsistent lengths",
            ));
        }
        Ok(self.distance.len())
    }

    fn with_capacity(rows: usize) -> Result<Self, String> {
        let Some(index_capacity) = rows.checked_mul(2) else {
            return Err("neighbor-table output exceeds addressable memory".to_owned());
        };
        let mut indices = Vec::new();
        let mut distance = Vec::new();
        indices
            .try_reserve_exact(index_capacity)
            .map_err(|_| "neighbor-table index allocation failed".to_owned())?;
        distance
            .try_reserve_exact(rows)
            .map_err(|_| "neighbor-table distance allocation failed".to_owned())?;
        Ok(Self { indices, distance })
    }

    pub(crate) fn from_pairs(pairs: Vec<pdbiox::NeighborPair>) -> Result<Self, String> {
        let mut table = Self::with_capacity(pairs.len())?;
        for pair in pairs {
            table.indices.extend([pair.first, pair.second]);
            table.distance.push(pair.distance_squared.sqrt());
        }
        Ok(table)
    }

    fn from_nearest(query: u32, neighbors: Vec<(u32, f32)>) -> Result<Self, String> {
        let mut table = Self::with_capacity(neighbors.len())?;
        for (target, distance_squared) in neighbors {
            table.indices.extend([query, target]);
            table.distance.push(distance_squared.sqrt());
        }
        Ok(table)
    }
}

enum SearchProfile {
    Backend(pdbiox::SpatialBackend),
    Options(pdbiox::SpatialSearchOptions),
}

pub(crate) enum SpatialBindingError {
    Input(String),
    Spatial(pdbiox::SpatialError),
    Allocation(String),
}

impl SpatialBindingError {
    pub(crate) fn into_pyerr(self) -> PyErr {
        match self {
            Self::Input(message) => PyValueError::new_err(message),
            Self::Spatial(error) => crate::spatial_index::spatial_error(error),
            Self::Allocation(message) => PyMemoryError::new_err(message),
        }
    }
}

/// Finds all unique neighbour pairs within `cutoff`.
#[pyfunction]
#[pyo3(signature = (positions, cutoff, *, left=None, right=None, backend=PySpatialBackend::Auto, cell=None))]
fn neighbor_pairs(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    cutoff: f32,
    left: Option<PyReadonlyArray1<'_, u32>>,
    right: Option<PyReadonlyArray1<'_, u32>>,
    backend: PySpatialBackend,
    cell: Option<&PyUnitCell>,
) -> PyResult<PyNeighborTable> {
    let positions = borrowed_coordinates(&positions)?;
    let left = index_slice(left.as_ref())?;
    let right = index_slice(right.as_ref())?;
    let periodic = periodic_box(cell)?;
    py.detach(move || {
        run_neighbor_pairs(
            positions,
            cutoff,
            left,
            right,
            SearchProfile::Backend(backend.into()),
            periodic.as_ref(),
        )
    })
    .map_err(SpatialBindingError::into_pyerr)
}

/// Finds all unique neighbour pairs with an explicit native search profile.
#[pyfunction]
#[pyo3(signature = (positions, cutoff, options, *, left=None, right=None, cell=None))]
fn neighbor_pairs_with_options(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    cutoff: f32,
    options: PySpatialSearchOptions,
    left: Option<PyReadonlyArray1<'_, u32>>,
    right: Option<PyReadonlyArray1<'_, u32>>,
    cell: Option<&PyUnitCell>,
) -> PyResult<PyNeighborTable> {
    let positions = borrowed_coordinates(&positions)?;
    let left = index_slice(left.as_ref())?;
    let right = index_slice(right.as_ref())?;
    let periodic = periodic_box(cell)?;
    py.detach(move || {
        run_neighbor_pairs(
            positions,
            cutoff,
            left,
            right,
            SearchProfile::Options(options.inner()),
            periodic.as_ref(),
        )
    })
    .map_err(SpatialBindingError::into_pyerr)
}

/// Selects query atom indices within `cutoff` of target atom indices.
#[pyfunction]
#[pyo3(signature = (positions, targets, cutoff, *, query=None, backend=PySpatialBackend::Auto, cell=None))]
fn atoms_within<'py>(
    py: Python<'py>,
    positions: PyReadonlyArray2<'_, f32>,
    targets: PyReadonlyArray1<'_, u32>,
    cutoff: f32,
    query: Option<PyReadonlyArray1<'_, u32>>,
    backend: PySpatialBackend,
    cell: Option<&PyUnitCell>,
) -> PyResult<Bound<'py, PyArray1<u32>>> {
    let positions = borrowed_coordinates(&positions)?;
    let targets = required_index_slice(&targets)?;
    let query = index_slice(query.as_ref())?;
    let periodic = periodic_box(cell)?;
    let indices = py
        .detach(move || {
            run_atoms_within(
                positions,
                targets,
                cutoff,
                query,
                SearchProfile::Backend(backend.into()),
                periodic.as_ref(),
            )
        })
        .map_err(SpatialBindingError::into_pyerr)?;
    Ok(readonly_index_result(py, indices))
}

/// Selects atom indices with an explicit native search profile.
#[pyfunction]
#[pyo3(signature = (positions, targets, cutoff, options, *, query=None, cell=None))]
fn atoms_within_with_options<'py>(
    py: Python<'py>,
    positions: PyReadonlyArray2<'_, f32>,
    targets: PyReadonlyArray1<'_, u32>,
    cutoff: f32,
    options: PySpatialSearchOptions,
    query: Option<PyReadonlyArray1<'_, u32>>,
    cell: Option<&PyUnitCell>,
) -> PyResult<Bound<'py, PyArray1<u32>>> {
    let positions = borrowed_coordinates(&positions)?;
    let targets = required_index_slice(&targets)?;
    let query = index_slice(query.as_ref())?;
    let periodic = periodic_box(cell)?;
    let indices = py
        .detach(move || {
            run_atoms_within(
                positions,
                targets,
                cutoff,
                query,
                SearchProfile::Options(options.inner()),
                periodic.as_ref(),
            )
        })
        .map_err(SpatialBindingError::into_pyerr)?;
    Ok(readonly_index_result(py, indices))
}

/// Finds the nearest target atoms for one query atom.
#[pyfunction]
#[pyo3(signature = (positions, targets, query, count, *, cell=None, periodic_options=None))]
fn nearest_neighbors(
    py: Python<'_>,
    positions: PyReadonlyArray2<'_, f32>,
    targets: PyReadonlyArray1<'_, u32>,
    query: u32,
    count: usize,
    cell: Option<&PyUnitCell>,
    periodic_options: Option<PyKdPeriodicOptions>,
) -> PyResult<PyNeighborTable> {
    let positions = borrowed_coordinates(&positions)?;
    let targets = required_index_slice(&targets)?;
    let periodic = periodic_box(cell)?;
    let periodic_options = match periodic_options {
        Some(options) => options.inner(),
        None => pdbiox::KdPeriodicOptions::default(),
    };
    py.detach(move || {
        validate_sorted_indices(targets, positions.len()).map_err(SpatialBindingError::Input)?;
        let neighbors = pdbiox::KdTree::build_with_options(
            positions,
            targets,
            periodic.as_ref(),
            periodic_options,
        )
        .and_then(|tree| tree.k_nearest(query, count))
        .map_err(SpatialBindingError::Spatial)?;
        PyNeighborTable::from_nearest(query, neighbors).map_err(SpatialBindingError::Allocation)
    })
    .map_err(SpatialBindingError::into_pyerr)
}

fn run_neighbor_pairs(
    positions: &[[f32; 3]],
    cutoff: f32,
    left: Option<&[u32]>,
    right: Option<&[u32]>,
    profile: SearchProfile,
    periodic: Option<&pdbiox::PeriodicBox>,
) -> Result<PyNeighborTable, SpatialBindingError> {
    let left = array_selection(left, positions.len()).map_err(SpatialBindingError::Input)?;
    let right = array_selection(right, positions.len()).map_err(SpatialBindingError::Input)?;
    let pairs = match profile {
        SearchProfile::Backend(backend) => pdbiox::pairs_within(
            positions,
            &left,
            &right,
            cutoff,
            backend,
            periodic,
            &crate::core::execution::default_context(),
        ),
        SearchProfile::Options(options) => pdbiox::pairs_within_with_options(
            positions,
            &left,
            &right,
            cutoff,
            options,
            periodic,
            &crate::core::execution::default_context(),
        ),
    }
    .map_err(SpatialBindingError::Spatial)?;
    PyNeighborTable::from_pairs(pairs).map_err(SpatialBindingError::Allocation)
}

fn run_atoms_within(
    positions: &[[f32; 3]],
    targets: &[u32],
    cutoff: f32,
    query: Option<&[u32]>,
    profile: SearchProfile,
    periodic: Option<&pdbiox::PeriodicBox>,
) -> Result<Vec<u32>, SpatialBindingError> {
    let query = array_selection(query, positions.len()).map_err(SpatialBindingError::Input)?;
    let targets =
        array_selection(Some(targets), positions.len()).map_err(SpatialBindingError::Input)?;
    let selected = match profile {
        SearchProfile::Backend(backend) => pdbiox::within(
            positions,
            &query,
            &targets,
            cutoff,
            backend,
            periodic,
            &crate::core::execution::default_context(),
        ),
        SearchProfile::Options(options) => pdbiox::within_with_options(
            positions,
            &query,
            &targets,
            cutoff,
            options,
            periodic,
            &crate::core::execution::default_context(),
        ),
    }
    .map_err(SpatialBindingError::Spatial)?;
    Ok(selected.into_iter().collect())
}

fn index_slice<'a>(indices: Option<&'a PyReadonlyArray1<'_, u32>>) -> PyResult<Option<&'a [u32]>> {
    indices
        .map(|array| {
            array.as_slice().map_err(|_| {
                PyValueError::new_err(
                    "atom indices must be C-contiguous uint32 arrays; pass a contiguous array explicitly",
                )
            })
        })
        .transpose()
}

pub(crate) fn required_index_slice<'a>(
    indices: &'a PyReadonlyArray1<'_, u32>,
) -> PyResult<&'a [u32]> {
    indices.as_slice().map_err(|_| {
        PyValueError::new_err(
            "atom indices must be C-contiguous uint32 arrays; pass a contiguous array explicitly",
        )
    })
}

fn array_selection(indices: Option<&[u32]>, count: usize) -> Result<pdbiox::AtomSelection, String> {
    let Some(indices) = indices else {
        return normalize_selection(None, count);
    };
    validate_sorted_indices(indices, count)?;
    Ok(pdbiox::AtomSelection::from_sorted(indices.to_vec()))
}

pub(crate) fn validate_sorted_indices(indices: &[u32], count: usize) -> Result<(), String> {
    let mut previous = None;
    for &index in indices {
        if index as usize >= count {
            return Err("atom index is outside the coordinate array".to_owned());
        }
        if let Some(previous) = previous
            && index <= previous
        {
            return Err("atom index arrays must be sorted and unique".to_owned());
        }
        previous = Some(index);
    }
    Ok(())
}

fn periodic_box(cell: Option<&PyUnitCell>) -> PyResult<Option<pdbiox::PeriodicBox>> {
    cell.map(|value| pdbiox::PeriodicBox::from_cell(value.cell))
        .transpose()
        .map_err(super::value_error)
}

fn readonly_indices<'py>(
    owner: &Bound<'py, PyNeighborTable>,
    rows: usize,
    pointer: *const u32,
) -> Bound<'py, PyArray2<u32>> {
    // SAFETY: `pointer` addresses `owner.indices`, whose row count was checked
    // against the exact two-column layout. The Python base retains `owner`.
    let view = unsafe { ArrayView2::from_shape_ptr((rows, 2), pointer) };
    // SAFETY: `owner` retains the private contiguous allocation for this view's
    // lifetime, and the frozen class exposes no Rust mutation of that buffer.
    let result = unsafe { PyArray2::borrow_from_array(&view, owner.clone().into_any()) };
    let _readonly = result.readwrite().make_nonwriteable();
    result
}

fn readonly_distances<'py>(
    owner: &Bound<'py, PyNeighborTable>,
    pointer: *const f32,
) -> Bound<'py, PyArray1<f32>> {
    let rows = owner.borrow().distance.len();
    // SAFETY: `pointer` addresses the private contiguous distance column, and
    // `owner` becomes the NumPy base object retaining it for the view lifetime.
    let view = unsafe { ArrayView1::from_shape_ptr(rows, pointer) };
    let result = unsafe { PyArray1::borrow_from_array(&view, owner.clone().into_any()) };
    let _readonly = result.readwrite().make_nonwriteable();
    result
}

fn readonly_index_result(py: Python<'_>, indices: Vec<u32>) -> Bound<'_, PyArray1<u32>> {
    let result = Array1::from_vec(indices).into_pyarray(py);
    let _readonly = result.readwrite().make_nonwriteable();
    result
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyNeighborTable>()?;
    module.add_function(wrap_pyfunction!(neighbor_pairs, module)?)?;
    module.add_function(wrap_pyfunction!(neighbor_pairs_with_options, module)?)?;
    module.add_function(wrap_pyfunction!(atoms_within, module)?)?;
    module.add_function(wrap_pyfunction!(atoms_within_with_options, module)?)?;
    module.add_function(wrap_pyfunction!(nearest_neighbors, module)?)?;
    Ok(())
}

#[cfg(test)]
#[path = "spatial_arrays_tests.rs"]
mod tests;
