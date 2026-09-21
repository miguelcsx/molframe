//! Python ownership and operations for contact tables.

use super::{PySelection, PyStructure};
use numpy::ndarray::ArrayView1;
use numpy::{PyArray1, PyArrayMethods};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyCapsule};
use std::sync::Arc;

#[derive(Debug)]
#[pyclass(name = "ContactTable", frozen, skip_from_py_object)]
pub(crate) struct PyContactTable {
    len: usize,
    contacts: Arc<molframe::analysis::ContactTable>,
    first: Py<PyArray1<u32>>,
    second: Py<PyArray1<u32>>,
    distance: Py<PyArray1<f32>>,
}

#[derive(Clone, Debug)]
#[pyclass(frozen, skip_from_py_object)]
struct ContactTableOwner {
    table: Arc<molframe::analysis::ContactTable>,
}

impl PyContactTable {
    pub(crate) fn from_native(
        py: Python<'_>,
        rows: molframe::analysis::ContactTable,
    ) -> PyResult<Self> {
        Self::from_shared(py, Arc::new(rows))
    }

    pub(crate) fn from_shared(
        py: Python<'_>,
        rows: Arc<molframe::analysis::ContactTable>,
    ) -> PyResult<Self> {
        let contacts = rows.clone();
        let owner = Bound::new(py, ContactTableOwner { table: rows })?;
        let (len, first, second, distance) = {
            let borrowed = owner.borrow();
            (
                borrowed.table.len(),
                borrowed.table.first().as_ptr().cast::<u32>(),
                borrowed.table.second().as_ptr().cast::<u32>(),
                borrowed.table.distances().as_ptr(),
            )
        };
        // SAFETY: typed indices are transparent `u32` values. The immutable
        // owner is installed as every NumPy base object and retains all columns.
        let first = unsafe { ArrayView1::from_shape_ptr(len, first) };
        // SAFETY: identical ownership and layout argument as `first`.
        let second = unsafe { ArrayView1::from_shape_ptr(len, second) };
        // SAFETY: `distance` points into the same retained immutable table.
        let distance = unsafe { ArrayView1::from_shape_ptr(len, distance) };
        // SAFETY: each borrowed array receives a strong reference to `owner`.
        let first = unsafe { PyArray1::borrow_from_array(&first, owner.clone().into_any()) };
        // SAFETY: each borrowed array receives a strong reference to `owner`.
        let second = unsafe { PyArray1::borrow_from_array(&second, owner.clone().into_any()) };
        // SAFETY: each borrowed array receives a strong reference to `owner`.
        let distance = unsafe { PyArray1::borrow_from_array(&distance, owner.into_any()) };
        first.readwrite().make_nonwriteable();
        second.readwrite().make_nonwriteable();
        distance.readwrite().make_nonwriteable();
        Ok(Self {
            len,
            contacts,
            first: first.unbind(),
            second: second.unbind(),
            distance: distance.unbind(),
        })
    }
}

#[pymethods]
impl PyContactTable {
    fn __len__(&self) -> usize {
        self.len
    }

    #[getter]
    fn first(&self, py: Python<'_>) -> Py<PyArray1<u32>> {
        self.first.clone_ref(py)
    }

    #[getter]
    fn second(&self, py: Python<'_>) -> Py<PyArray1<u32>> {
        self.second.clone_ref(py)
    }

    #[getter]
    fn distance(&self, py: Python<'_>) -> Py<PyArray1<f32>> {
        self.distance.clone_ref(py)
    }

    #[pyo3(signature = (_requested_schema=None))]
    fn __arrow_c_stream__<'py>(
        &self,
        py: Python<'py>,
        _requested_schema: Option<&Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyCapsule>> {
        let stream = molframe::interop::ContactArrowTable::new(self.contacts.clone())
            .arrow_stream()
            .map_err(|error| pyo3::exceptions::PyRuntimeError::new_err(error.to_string()))?;
        PyCapsule::new_with_value(py, stream.into_ffi(), c"arrow_array_stream")
    }
}

#[pyfunction]
#[pyo3(signature = (value, cutoff, *, backend="auto"))]
pub(crate) fn atom_contacts(
    py: Python<'_>,
    value: &Bound<'_, PyAny>,
    cutoff: f32,
    backend: &str,
) -> PyResult<PyContactTable> {
    let backend = parse_backend(backend)?;
    let rows = if let Ok(structure) = value.extract::<PyRef<'_, PyStructure>>() {
        let structure = structure.inner.clone();
        py.detach(move || {
            molframe::analysis::atom_contacts(
                structure.engine(),
                cutoff,
                backend,
                &molframe::ExecutionContext::default(),
            )
        })
    } else if let Ok(selection) = value.extract::<PyRef<'_, PySelection>>() {
        let structure = selection.parent.inner.clone();
        let indices = selection.indices.clone();
        py.detach(move || {
            let selected = molframe::engine::core::AtomSelection::from_sorted(indices);
            molframe::analysis::atom_contacts_between(
                structure.engine(),
                &selected,
                &selected,
                cutoff,
                backend,
                &molframe::ExecutionContext::default(),
            )
        })
    } else {
        return Err(pyo3::exceptions::PyTypeError::new_err(
            "value must be a Structure or structure-bound Selection",
        ));
    }
    .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))?;
    PyContactTable::from_native(py, rows)
}

fn parse_backend(value: &str) -> PyResult<molframe::spatial::SpatialBackend> {
    match value {
        "auto" => Ok(molframe::spatial::SpatialBackend::Auto),
        "cell" => Ok(molframe::spatial::SpatialBackend::CellList),
        "kd_tree" => Ok(molframe::spatial::SpatialBackend::KdTree),
        "brute_force" => Ok(molframe::spatial::SpatialBackend::BruteForce),
        _ => Err(pyo3::exceptions::PyValueError::new_err(
            "backend must be 'auto', 'cell', 'kd_tree', or 'brute_force'",
        )),
    }
}
