//! Central facade handles and stateless namespace functions.

use crate::hierarchy::{PyAtoms, PyChains, PyModels, PyResidues};
#[cfg(feature = "geometry")]
use numpy::PyReadonlyArray2;
#[cfg(feature = "geometry")]
use numpy::PyUntypedArrayMethods;
#[cfg(feature = "analysis")]
use numpy::ndarray::ArrayView1;
use numpy::ndarray::ArrayView2;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods};
use pyo3::prelude::*;
use pyo3::pybacked::PyBackedBytes;
use pyo3::types::PyAny;
#[cfg(feature = "analysis")]
use pyo3::types::PyCapsule;
use std::fmt;
use std::path::PathBuf;
#[cfg(feature = "analysis")]
use std::sync::Arc;

#[derive(Clone, Debug)]
#[pyclass(name = "Structure", frozen, skip_from_py_object)]
pub(crate) struct PyStructure {
    pub(crate) inner: molframe::Structure,
}
#[derive(Clone, Debug)]
#[pyclass(frozen, skip_from_py_object)]
struct SnapshotOwner {
    structure: molframe::Structure,
}
impl PyStructure {
    pub(crate) const fn new(inner: molframe::Structure) -> Self {
        Self { inner }
    }

    fn coordinates_view<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let owner = Bound::new(
            py,
            SnapshotOwner {
                structure: self.inner.clone(),
            },
        )?;
        let (rows, pointer) = {
            let snapshot = owner.borrow();
            let coordinates = snapshot.structure.coordinates();
            (coordinates.len(), coordinates.as_ptr().cast::<f32>())
        };
        rows.checked_mul(3).ok_or_else(|| {
            pyo3::exceptions::PyOverflowError::new_err("coordinate array shape overflows usize")
        })?;
        // SAFETY: `SnapshotOwner` retains the immutable Arc-backed snapshot,
        // `[f32; 3]` is contiguous, and the resulting array is made read-only.
        let view = unsafe { ArrayView2::from_shape_ptr((rows, 3), pointer) };
        // SAFETY: the NumPy base object is exactly the owner used to obtain the
        // pointer, so Python cannot outlive or mutate the backing allocation.
        let array = unsafe { PyArray2::borrow_from_array(&view, owner.into_any()) };
        array.readwrite().make_nonwriteable();
        Ok(array)
    }
}
#[pymethods]
impl PyStructure {
    fn __repr__(&self) -> String {
        self.inner.to_string()
    }

    #[getter]
    fn coordinates<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        self.coordinates_view(py)
    }

    #[getter]
    fn atoms(&self) -> PyAtoms {
        PyAtoms::all(self.clone())
    }

    #[getter]
    fn residues(&self) -> PyResidues {
        PyResidues::all(self.clone())
    }

    #[getter]
    fn chains(&self) -> PyChains {
        PyChains::all(self.clone())
    }

    #[getter]
    fn models(&self) -> PyModels {
        PyModels::new(self.clone())
    }

    #[pyo3(signature = (query, *, policy=None))]
    fn select(&self, query: &str, policy: Option<&Bound<'_, PyAny>>) -> PyResult<PySelection> {
        if policy.is_some() {
            return Err(pyo3::exceptions::PyNotImplementedError::new_err(
                "custom policy objects are not yet accepted by this binding",
            ));
        }
        selection(self, query)
    }

    fn edit(&self) -> PyStructureEditor {
        PyStructureEditor {
            inner: Some(self.inner.edit()),
        }
    }
}
#[derive(Debug)]
#[pyclass(name = "StructureEditor", skip_from_py_object)]
pub(crate) struct PyStructureEditor {
    inner: Option<molframe::StructureEditor>,
}
#[pymethods]
impl PyStructureEditor {
    fn rename_chain(&mut self, chain: u32, label: &str) -> PyResult<()> {
        let Some(editor) = self.inner.as_mut() else {
            return Err(pyo3::exceptions::PyRuntimeError::new_err(
                "editor has already been finished",
            ));
        };
        editor
            .rename_chain(molframe::ChainIndex::new(chain), label)
            .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))
    }

    fn finish(&mut self) -> PyResult<PyStructure> {
        let Some(editor) = self.inner.take() else {
            return Err(pyo3::exceptions::PyRuntimeError::new_err(
                "editor has already been finished",
            ));
        };
        editor
            .finish()
            .map(PyStructure::new)
            .map_err(|diagnostics| {
                pyo3::exceptions::PyValueError::new_err(format!("{diagnostics:?}"))
            })
    }
}
#[derive(Clone, Debug)]
#[pyclass(name = "Selection", frozen, skip_from_py_object)]
pub(crate) struct PySelection {
    parent: PyStructure,
    selection: molframe::Selection,
    indices: Vec<u32>,
}
#[pymethods]
impl PySelection {
    fn __len__(&self) -> usize {
        match usize::try_from(self.selection.len()) {
            Ok(value) => value,
            Err(_) => usize::MAX,
        }
    }

    fn is_stale_for(&self, structure: &PyStructure) -> bool {
        self.selection.is_stale_for(&structure.inner)
    }

    #[getter]
    fn indices<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u32>> {
        self.indices.clone().into_pyarray(py)
    }

    #[pyo3(signature = (*, model=0))]
    fn to_coordinates<'py>(
        &self,
        py: Python<'py>,
        model: u32,
    ) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let coordinates = self
            .selection
            .to_coordinates(molframe::ModelIndex::new(model));
        let rows = coordinates.len();
        let flat = coordinates.into_iter().flatten().collect::<Vec<_>>();
        let array = flat.into_pyarray(py);
        array.reshape((rows, 3))
    }

    fn __or__(&self, other: &Self) -> Self {
        let combined = &self.selection | &other.selection;
        Self::from_native(self.parent.clone(), combined)
    }

    fn __and__(&self, other: &Self) -> Self {
        let combined = &self.selection & &other.selection;
        Self::from_native(self.parent.clone(), combined)
    }

    fn __sub__(&self, other: &Self) -> Self {
        let combined = &self.selection - &other.selection;
        Self::from_native(self.parent.clone(), combined)
    }
}

impl PySelection {
    fn from_native(parent: PyStructure, selection: molframe::Selection) -> Self {
        let indices = selection.atoms().map(|atom| atom.index().get()).collect();
        Self {
            parent,
            selection,
            indices,
        }
    }
}

#[derive(Clone, Debug)]
#[pyclass(name = "Query", frozen, skip_from_py_object)]
pub(crate) struct PyQuery {
    compiled: molframe::Query,
}

impl PyQuery {
    pub(crate) fn from_native(compiled: molframe::Query) -> Self {
        Self { compiled }
    }

    pub(crate) fn native(&self) -> &molframe::Query {
        &self.compiled
    }
}

#[pymethods]
impl PyQuery {
    #[new]
    fn new(source: &str) -> PyResult<Self> {
        let compiled = molframe::Query::compile(source)
            .map_err(|findings| findings_error(molframe::Findings::from(findings)))?;
        Ok(Self { compiled })
    }

    fn select(&self, structure: &PyStructure) -> PyResult<PySelection> {
        use molframe::QueryStructure;

        let evaluation = structure
            .inner
            .select_query(&self.compiled, &molframe::AnalysisPolicy::default())
            .map_err(findings_error)?;
        let view = structure.inner.engine().view_of(evaluation.selection);
        Ok(PySelection::from_native(
            structure.clone(),
            molframe::Selection::from(view),
        ))
    }

    fn __and__(&self, other: &Self) -> Self {
        Self::from_native(self.compiled.clone() & other.compiled.clone())
    }

    fn __or__(&self, other: &Self) -> Self {
        Self::from_native(self.compiled.clone() | other.compiled.clone())
    }

    fn __invert__(&self) -> Self {
        Self::from_native(!self.compiled.clone())
    }
}

struct PythonBytes(PyBackedBytes);

impl AsRef<[u8]> for PythonBytes {
    fn as_ref(&self) -> &[u8] {
        self.0.as_ref()
    }
}

impl fmt::Debug for PythonBytes {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("PythonBytes")
            .field("len", &self.0.len())
            .finish()
    }
}

#[derive(Clone, Debug)]
#[pyclass(name = "Reader", frozen, skip_from_py_object)]
pub(crate) struct PyReader {
    input: molframe::InputBuffer,
    name: Option<Box<str>>,
}

#[pymethods]
impl PyReader {
    #[new]
    #[pyo3(signature = (data, *, name=None))]
    fn new(data: PyBackedBytes, name: Option<&str>) -> Self {
        Self {
            input: molframe::InputBuffer::from_owner(PythonBytes(data)),
            name: name.map(Into::into),
        }
    }

    fn read(&self, py: Python<'_>) -> PyResult<PyStructure> {
        let input = self.input.clone();
        let name = self.name.clone();
        py.detach(move || {
            molframe::read_buffer(&input, name.as_deref(), &molframe::ReadOptions::new())
        })
        .map(|(structure, _)| PyStructure::new(structure))
        .map_err(findings_error)
    }

    #[getter]
    fn byte_length(&self) -> usize {
        self.input.len()
    }
}

#[pyfunction]
#[pyo3(signature = (source, *, name=None))]
pub(crate) fn read(
    py: Python<'_>,
    source: &Bound<'_, PyAny>,
    name: Option<&str>,
) -> PyResult<PyStructure> {
    if let Ok(path) = source.extract::<PathBuf>() {
        return py
            .detach(move || molframe::read(path))
            .map(PyStructure::new)
            .map_err(findings_error);
    }
    let bytes = source.extract::<PyBackedBytes>()?;
    PyReader::new(bytes, name).read(py)
}

fn selection(structure: &PyStructure, source: &str) -> PyResult<PySelection> {
    structure
        .inner
        .select(source, &molframe::AnalysisPolicy::default())
        .map(|selection| PySelection::from_native(structure.clone(), selection))
        .map_err(findings_error)
}

#[cfg(feature = "geometry")]
pub(crate) fn coordinates<'a>(array: &'a PyReadonlyArray2<'_, f32>) -> PyResult<&'a [[f32; 3]]> {
    let shape = array.shape();
    if shape.len() != 2 || shape[1] != 3 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "coordinates must have shape (n, 3)",
        ));
    }
    let contiguous = array.as_slice().map_err(|_| {
        pyo3::exceptions::PyValueError::new_err(
            "coordinates must be C-contiguous; call numpy.ascontiguousarray",
        )
    })?;
    Ok(contiguous.as_chunks::<3>().0)
}

#[pyfunction]
#[cfg(feature = "geometry")]
pub(crate) fn centroid(array: PyReadonlyArray2<'_, f32>) -> PyResult<Option<[f64; 3]>> {
    Ok(molframe::geometry::centroid(coordinates(&array)?))
}

#[pyfunction]
#[cfg(feature = "geometry")]
pub(crate) fn rmsd(
    mobile: PyReadonlyArray2<'_, f32>,
    reference: PyReadonlyArray2<'_, f32>,
) -> PyResult<f64> {
    molframe::geometry::rmsd(coordinates(&mobile)?, coordinates(&reference)?)
        .map_err(|error| pyo3::exceptions::PyValueError::new_err(format!("{error:?}")))
}

#[pyfunction]
#[cfg(feature = "geometry")]
pub(crate) fn distance_matrix<'py>(
    py: Python<'py>,
    array: PyReadonlyArray2<'_, f32>,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let matrix = molframe::geometry::distance_matrix(coordinates(&array)?)
        .map_err(|error| pyo3::exceptions::PyMemoryError::new_err(error.to_string()))?;
    let rows = matrix.rows();
    matrix.into_values().into_pyarray(py).reshape((rows, rows))
}

#[derive(Debug)]
#[cfg(feature = "analysis")]
#[pyclass(name = "ContactTable", frozen, skip_from_py_object)]
pub(crate) struct PyContactTable {
    len: usize,
    contacts: Arc<molframe::analysis::ContactTable>,
    first: Py<PyArray1<u32>>,
    second: Py<PyArray1<u32>>,
    distance: Py<PyArray1<f32>>,
}

#[derive(Clone, Debug)]
#[cfg(feature = "analysis")]
#[pyclass(frozen, skip_from_py_object)]
struct ContactTableOwner {
    table: Arc<molframe::analysis::ContactTable>,
}

#[cfg(feature = "analysis")]
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
#[cfg(feature = "analysis")]
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
#[cfg(feature = "analysis")]
#[pyo3(signature = (value, cutoff, *, backend="auto"))]
pub(crate) fn atom_contacts(
    py: Python<'_>,
    value: &Bound<'_, PyAny>,
    cutoff: f32,
    backend: &str,
) -> PyResult<PyContactTable> {
    let backend = match backend {
        "auto" => molframe::spatial::SpatialBackend::Auto,
        "cell" => molframe::spatial::SpatialBackend::CellList,
        "kd_tree" => molframe::spatial::SpatialBackend::KdTree,
        "brute_force" => molframe::spatial::SpatialBackend::BruteForce,
        _ => {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "backend must be 'auto', 'cell', 'kd_tree', or 'brute_force'",
            ));
        }
    };
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

fn findings_error(findings: molframe::Findings) -> PyErr {
    pyo3::exceptions::PyValueError::new_err(findings.to_string())
}
