//! Shared payload wrappers and zero-copy `NumPy` projections.

use super::metadata::{PyChunkDescriptor, index_error, value_error};
use crate::bonds::{PyBondOrder, PyBondProvenance};
use crate::chemistry::PyElement;
use crate::core_columns::PyPresence;
use crate::core_storage::PyCoordinateBlock;
use crate::core_values::PySymbolId;
use crate::index::PyLocalRow;
use crate::index::{PyDatasetId, PyLogicalRow};
use crate::structure::PyStructure;
use numpy::ndarray::{ArrayView1, ArrayView2};
use numpy::{PyArray1, PyArray2, PyArrayMethods};
use pyo3::PyClass;
use pyo3::prelude::*;

#[pyclass(name = "PropertyKind", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyPropertyKind {
    Boolean,
    Integer,
    Real,
    Symbol,
}

impl From<molframe::PropertyKind> for PyPropertyKind {
    fn from(value: molframe::PropertyKind) -> Self {
        match value {
            molframe::PropertyKind::Boolean => Self::Boolean,
            molframe::PropertyKind::Integer => Self::Integer,
            molframe::PropertyKind::Real => Self::Real,
            molframe::PropertyKind::Symbol => Self::Symbol,
        }
    }
}

#[pyclass(name = "PropertyValue", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyPropertyValue(molframe::PropertyValue);

#[pyclass(name = "AtomEndpoint", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyAtomEndpoint(molframe::AtomEndpoint);

#[pymethods]
impl PyAtomEndpoint {
    #[new]
    fn new(dataset: PyDatasetId, row: PyLogicalRow) -> Self {
        Self(molframe::AtomEndpoint::new(dataset.0, row.0))
    }

    #[getter]
    fn dataset(&self) -> PyDatasetId {
        PyDatasetId(self.0.dataset())
    }

    #[getter]
    fn row(&self) -> PyLogicalRow {
        PyLogicalRow(self.0.row())
    }
}

#[pyclass(name = "BondChunkRecord", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyBondChunkRecord(molframe::BondChunkRecord);

#[pymethods]
impl PyBondChunkRecord {
    #[getter]
    fn atom_a(&self) -> PyAtomEndpoint {
        PyAtomEndpoint(self.0.atom_a)
    }

    #[getter]
    fn atom_b(&self) -> PyAtomEndpoint {
        PyAtomEndpoint(self.0.atom_b)
    }

    #[getter]
    fn order(&self) -> PyBondOrder {
        self.0.order.into()
    }

    #[getter]
    fn provenance(&self) -> PyBondProvenance {
        self.0.provenance.into()
    }
}

#[pyclass(name = "BondChunk", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyBondChunk(pub(crate) molframe::BondChunk);

#[pymethods]
impl PyBondChunk {
    #[getter]
    fn descriptor(&self) -> PyChunkDescriptor {
        PyChunkDescriptor(self.0.descriptor())
    }

    #[getter]
    fn atom_dataset(&self) -> PyDatasetId {
        PyDatasetId(self.0.atom_dataset())
    }

    #[getter]
    fn structure(&self) -> PyStructure {
        PyStructure::new(self.0.structure().clone())
    }

    fn record(&self, local: PyLocalRow) -> PyResult<PyBondChunkRecord> {
        self.0
            .record(local.0)
            .map(PyBondChunkRecord)
            .map_err(index_error)
    }

    fn minimum_host_bytes(&self) -> PyResult<u64> {
        self.0.minimum_host_bytes().map_err(value_error)
    }
}

#[pymethods]
impl PyPropertyValue {
    #[getter]
    fn kind(&self) -> PyPropertyKind {
        match self.0 {
            molframe::PropertyValue::Boolean(_, _) => PyPropertyKind::Boolean,
            molframe::PropertyValue::Integer(_, _) => PyPropertyKind::Integer,
            molframe::PropertyValue::Real(_, _) => PyPropertyKind::Real,
            molframe::PropertyValue::Symbol(_, _) => PyPropertyKind::Symbol,
        }
    }

    #[getter]
    fn presence(&self) -> PyPresence {
        match self.0 {
            molframe::PropertyValue::Boolean(_, presence)
            | molframe::PropertyValue::Integer(_, presence)
            | molframe::PropertyValue::Real(_, presence)
            | molframe::PropertyValue::Symbol(_, presence) => presence.into(),
        }
    }

    #[getter]
    fn boolean(&self) -> Option<bool> {
        match self.0 {
            molframe::PropertyValue::Boolean(value, _) => Some(value),
            _ => None,
        }
    }

    #[getter]
    fn integer(&self) -> Option<i64> {
        match self.0 {
            molframe::PropertyValue::Integer(value, _) => Some(value),
            _ => None,
        }
    }

    #[getter]
    fn real(&self) -> Option<f64> {
        match self.0 {
            molframe::PropertyValue::Real(value, _) => Some(value),
            _ => None,
        }
    }

    #[getter]
    fn symbol(&self) -> Option<PySymbolId> {
        match self.0 {
            molframe::PropertyValue::Symbol(value, _) => Some(PySymbolId(value)),
            _ => None,
        }
    }
}

#[pyclass(name = "StructureChunk", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyStructureChunk(pub(crate) molframe::StructureChunk);

#[pymethods]
impl PyStructureChunk {
    #[getter]
    fn descriptor(&self) -> PyChunkDescriptor {
        PyChunkDescriptor(self.0.descriptor())
    }

    #[getter]
    fn positions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        structure_positions(py, self.clone())
    }

    fn element(&self, local: PyLocalRow) -> Option<PyElement> {
        self.0.atoms().element(local.0.get()).map(PyElement)
    }

    fn atom_site_id(&self, local: PyLocalRow) -> Option<u32> {
        self.0.atoms().atom_site_id(local.0.get())
    }
}

#[pyclass(name = "PropertyChunk", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPropertyChunk(pub(crate) molframe::PropertyChunk);

#[pymethods]
impl PyPropertyChunk {
    #[getter]
    fn descriptor(&self) -> PyChunkDescriptor {
        PyChunkDescriptor(self.0.descriptor())
    }

    #[getter]
    fn name(&self) -> &str {
        self.0.name()
    }

    #[getter]
    fn kind(&self) -> Option<PyPropertyKind> {
        self.0.kind().map(Into::into)
    }

    fn value(&self, local: PyLocalRow) -> PyResult<PyPropertyValue> {
        self.0
            .value(local.0)
            .map(PyPropertyValue)
            .map_err(index_error)
    }

    fn integers<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyArray1<i64>>>> {
        shared_integers(py, self.clone())
    }

    fn reals<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyArray1<f64>>>> {
        shared_reals(py, self.clone())
    }
}

#[pyclass(name = "FrameChunk", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyFrameChunk(pub(crate) molframe::FrameChunk);

#[pymethods]
impl PyFrameChunk {
    #[new]
    fn new(
        descriptor: PyChunkDescriptor,
        coordinates: PyCoordinateBlock,
        start: u32,
        end: u32,
    ) -> PyResult<Self> {
        molframe::FrameChunk::shared(descriptor.0, coordinates.0, start..end)
            .map(Self)
            .map_err(value_error)
    }

    #[getter]
    fn descriptor(&self) -> PyChunkDescriptor {
        PyChunkDescriptor(self.0.descriptor())
    }

    #[getter]
    fn positions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        frame_positions(py, self.clone())
    }

    #[getter]
    fn coordinates(&self) -> PyCoordinateBlock {
        PyCoordinateBlock(self.0.coordinates().clone())
    }
}

fn structure_positions(
    py: Python<'_>,
    chunk: PyStructureChunk,
) -> PyResult<Bound<'_, PyArray2<f32>>> {
    let owner = Bound::new(py, chunk)?;
    let (rows, pointer) = {
        let borrowed = owner.borrow();
        let positions = borrowed.0.positions();
        (positions.len(), positions.as_ptr().cast::<f32>())
    };
    readonly_vectors(owner, rows, pointer)
}

fn frame_positions(py: Python<'_>, chunk: PyFrameChunk) -> PyResult<Bound<'_, PyArray2<f32>>> {
    let owner = Bound::new(py, chunk)?;
    let (rows, pointer) = {
        let borrowed = owner.borrow();
        let positions = borrowed.0.positions();
        (positions.len(), positions.as_ptr().cast::<f32>())
    };
    readonly_vectors(owner, rows, pointer)
}

fn readonly_vectors<T: PyClass>(
    owner: Bound<'_, T>,
    rows: usize,
    pointer: *const f32,
) -> PyResult<Bound<'_, PyArray2<f32>>> {
    rows.checked_mul(3)
        .ok_or_else(|| pyo3::exceptions::PyOverflowError::new_err("coordinate shape overflow"))?;
    // SAFETY: native chunks expose contiguous `[f32; 3]` rows and `owner`
    // retains their immutable backing allocation for the NumPy view lifetime.
    let view = unsafe { ArrayView2::from_shape_ptr((rows, 3), pointer) };
    // SAFETY: the borrowed array keeps `owner` as its Python base object.
    let array = unsafe { PyArray2::borrow_from_array(&view, owner.into_any()) };
    array.readwrite().make_nonwriteable();
    Ok(array)
}

fn shared_integers(
    py: Python<'_>,
    chunk: PyPropertyChunk,
) -> PyResult<Option<Bound<'_, PyArray1<i64>>>> {
    let owner = Bound::new(py, chunk)?;
    let Some((length, pointer)) = ({
        let borrowed = owner.borrow();
        borrowed
            .0
            .integers()
            .map(|values| (values.len(), values.as_ptr()))
    }) else {
        return Ok(None);
    };
    // SAFETY: `pointer` and `length` came from the immutable integer slice
    // retained by `owner`.
    let view = unsafe { ArrayView1::from_shape_ptr(length, pointer) };
    // SAFETY: the borrowed array keeps `owner` as its Python base object.
    let array = unsafe { PyArray1::borrow_from_array(&view, owner.into_any()) };
    array.readwrite().make_nonwriteable();
    Ok(Some(array))
}

fn shared_reals(
    py: Python<'_>,
    chunk: PyPropertyChunk,
) -> PyResult<Option<Bound<'_, PyArray1<f64>>>> {
    let owner = Bound::new(py, chunk)?;
    let Some((length, pointer)) = ({
        let borrowed = owner.borrow();
        borrowed
            .0
            .reals()
            .map(|values| (values.len(), values.as_ptr()))
    }) else {
        return Ok(None);
    };
    // SAFETY: `pointer` and `length` came from the immutable real slice
    // retained by `owner`.
    let view = unsafe { ArrayView1::from_shape_ptr(length, pointer) };
    // SAFETY: the borrowed array keeps `owner` as its Python base object.
    let array = unsafe { PyArray1::borrow_from_array(&view, owner.into_any()) };
    array.readwrite().make_nonwriteable();
    Ok(Some(array))
}
