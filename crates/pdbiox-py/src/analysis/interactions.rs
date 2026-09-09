//! Spatial contacts and secondary structure in coarse Rust calls.

use crate::contract::{PyAnalysis, analysis_with_value};
use crate::graph::PySpatialBackend;
use crate::query::PySelection;
use crate::structure::PyStructure;
use numpy::ndarray::ArrayView1;
use numpy::{PyArray1, PyArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "Contact", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyContact {
    first: u32,
    second: u32,
    distance: f32,
}

#[pymethods]
impl PyContact {
    #[new]
    fn new(first: u32, second: u32, distance: f32) -> Self {
        Self {
            first,
            second,
            distance,
        }
    }

    #[getter]
    const fn first(&self) -> u32 {
        self.first
    }
    #[getter]
    const fn second(&self) -> u32 {
        self.second
    }
    #[getter]
    const fn distance(&self) -> f32 {
        self.distance
    }

    fn __repr__(&self) -> String {
        format!(
            "Contact(first={}, second={}, distance={})",
            self.first, self.second, self.distance
        )
    }
}

/// Dense, columnar atom-contact results.
///
/// The three buffers are built together from the native contact kernel, so row
/// `i` always describes one contact across all columns. Their `NumPy` views use
/// this object as their base, retaining the Rust buffers without copying them.
#[pyclass(name = "ContactTable", frozen, skip_from_py_object)]
#[derive(Debug)]
pub(crate) struct PyContactTable {
    first: Vec<u32>,
    second: Vec<u32>,
    distance: Vec<f32>,
}

#[pymethods]
impl PyContactTable {
    #[getter]
    fn first<'py>(slf: &Bound<'py, Self>) -> Bound<'py, PyArray1<u32>> {
        let arrays = slf.borrow();
        readonly_array(slf, &arrays.first)
    }

    #[getter]
    fn second<'py>(slf: &Bound<'py, Self>) -> Bound<'py, PyArray1<u32>> {
        let arrays = slf.borrow();
        readonly_array(slf, &arrays.second)
    }

    #[getter]
    fn distance<'py>(slf: &Bound<'py, Self>) -> Bound<'py, PyArray1<f32>> {
        let arrays = slf.borrow();
        readonly_array(slf, &arrays.distance)
    }

    fn __len__(&self) -> usize {
        self.distance.len()
    }

    fn __repr__(&self) -> String {
        format!("ContactTable(rows={})", self.distance.len())
    }
}

#[pyclass(name = "SseKind", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PySseKind {
    AlphaHelix,
    Strand,
    Turn,
    Coil,
}

#[pyclass(name = "SecondaryStructure", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySecondaryStructure {
    residue: u32,
    kind: PySseKind,
}

/// Exact projection of the facade's `SseRecord` result.  The older
/// `SecondaryStructure` class is retained for the structure convenience method;
/// the free function below returns this lossless facade shape.
#[pyclass(name = "SseRecord", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PySseRecord {
    #[pyo3(get)]
    residue: u32,
    #[pyo3(get)]
    kind: PySseKind,
}

#[pyclass(name = "DsspOptions", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyDsspOptions(pdbiox::analysis::DsspOptions);

#[pymethods]
impl PyDsspOptions {
    #[new]
    fn new(
        electrostatic_prefactor: f64,
        hydrogen_bond_energy: f64,
        amide_hydrogen_distance: f32,
        minimum_sequence_separation: usize,
        helix_offset: usize,
        turn_offsets: (usize, usize),
    ) -> Self {
        Self(pdbiox::analysis::DsspOptions {
            electrostatic_prefactor,
            hydrogen_bond_energy,
            amide_hydrogen_distance,
            minimum_sequence_separation,
            helix_offset,
            turn_offsets: turn_offsets.0..=turn_offsets.1,
        })
    }
}

impl PyDsspOptions {
    pub(crate) fn native(&self) -> pdbiox::analysis::DsspOptions {
        self.0.clone()
    }

    pub(crate) fn from_native_parts(value: pdbiox::analysis::DsspOptions) -> Self {
        Self(value)
    }
}

#[pymethods]
impl PySecondaryStructure {
    #[getter]
    const fn residue(&self) -> u32 {
        self.residue
    }
    #[getter]
    const fn kind(&self) -> PySseKind {
        self.kind
    }
}

#[pymethods]
impl PyStructure {
    #[pyo3(signature = (cutoff, *, backend=PySpatialBackend::Auto))]
    fn contacts(
        &self,
        py: Python<'_>,
        cutoff: f32,
        backend: PySpatialBackend,
    ) -> PyResult<PyContactTable> {
        let structure = self.structure().clone();
        py.detach(move || {
            pdbiox::analysis::atom_contacts(
                &structure,
                cutoff,
                backend.into(),
                &crate::core::execution::default_context(),
            )
            .map(PyContactTable::from)
        })
        .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    fn secondary_structure(
        &self,
        py: Python<'_>,
        options: &PyDsspOptions,
    ) -> PyResult<Vec<PySecondaryStructure>> {
        let structure = self.structure().clone();
        let options = options.0.clone();
        py.detach(move || pdbiox::analysis::secondary_structure(&structure, &options))
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| PySecondaryStructure {
                        residue: value.residue.get(),
                        kind: value.kind.into(),
                    })
                    .collect()
            })
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }
}

/// Runs the facade atom-contact kernel without a Python-side atom loop.
#[pyfunction]
pub(crate) fn atom_contacts(
    py: Python<'_>,
    structure: &PyStructure,
    cutoff: f32,
    backend: PySpatialBackend,
) -> PyResult<PyContactTable> {
    let structure = structure.structure().clone();
    py.detach(move || {
        pdbiox::analysis::atom_contacts(
            &structure,
            cutoff,
            backend.into(),
            &crate::core::execution::default_context(),
        )
        .map(PyContactTable::from)
    })
    .map_err(value_error)
}

/// Runs the native contact kernel between two already compiled selections.
#[pyfunction]
pub(crate) fn atom_contacts_between(
    py: Python<'_>,
    structure: &PyStructure,
    left: &PySelection,
    right: &PySelection,
    cutoff: f32,
    backend: PySpatialBackend,
) -> PyResult<PyContactTable> {
    let structure = structure.structure().clone();
    let left = left.inner.clone();
    let right = right.inner.clone();
    py.detach(move || {
        pdbiox::analysis::atom_contacts_between(
            &structure,
            &left,
            &right,
            cutoff,
            backend.into(),
            &crate::core::execution::default_context(),
        )
        .map(PyContactTable::from)
    })
    .map_err(value_error)
}

/// Returns the facade's complete residue/state records for DSSP assignment.
#[pyfunction]
pub(crate) fn secondary_structure(
    py: Python<'_>,
    structure: &PyStructure,
    options: &PyDsspOptions,
) -> PyResult<Vec<PySseRecord>> {
    let structure = structure.structure().clone();
    let options = options.0.clone();
    py.detach(move || pdbiox::analysis::secondary_structure(&structure, &options))
        .map(|values| {
            values
                .into_iter()
                .map(|value| PySseRecord {
                    residue: value.residue.get(),
                    kind: value.kind.into(),
                })
                .collect()
        })
        .map_err(value_error)
}

impl From<pdbiox::analysis::Contact> for PyContact {
    fn from(value: pdbiox::analysis::Contact) -> Self {
        Self {
            first: value.first.get(),
            second: value.second.get(),
            distance: value.distance,
        }
    }
}

impl From<Vec<pdbiox::analysis::Contact>> for PyContactTable {
    fn from(contacts: Vec<pdbiox::analysis::Contact>) -> Self {
        let mut first = Vec::with_capacity(contacts.len());
        let mut second = Vec::with_capacity(contacts.len());
        let mut distance = Vec::with_capacity(contacts.len());
        for contact in contacts {
            first.push(contact.first.get());
            second.push(contact.second.get());
            distance.push(contact.distance);
        }
        Self {
            first,
            second,
            distance,
        }
    }
}

impl From<pdbiox::analysis::SseKind> for PySseKind {
    fn from(value: pdbiox::analysis::SseKind) -> Self {
        match value {
            pdbiox::analysis::SseKind::AlphaHelix => Self::AlphaHelix,
            pdbiox::analysis::SseKind::Strand => Self::Strand,
            pdbiox::analysis::SseKind::Turn => Self::Turn,
            pdbiox::analysis::SseKind::Coil => Self::Coil,
        }
    }
}

impl From<pdbiox::analysis::SseRecord> for PySseRecord {
    fn from(value: pdbiox::analysis::SseRecord) -> Self {
        Self {
            residue: value.residue.get(),
            kind: value.kind.into(),
        }
    }
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}

pub(crate) fn contact_analysis_to_py(
    py: Python<'_>,
    analysis: pdbiox::Analysis<Vec<pdbiox::analysis::Contact>>,
) -> PyResult<PyAnalysis> {
    let analysis = py.detach(move || analysis.map(PyContactTable::from));
    analysis_with_value(py, analysis, |py, table| Ok(Py::new(py, table)?.into_any()))
}

fn readonly_array<'py, T>(
    owner: &Bound<'py, PyContactTable>,
    values: &[T],
) -> Bound<'py, PyArray1<T>>
where
    T: numpy::Element,
{
    // SAFETY: `values` is one of `owner`'s private vectors. `owner` becomes
    // the NumPy base object, so it retains the contiguous allocation for the
    // full lifetime of the view; this frozen class exposes no Rust mutation.
    let view = unsafe { ArrayView1::from_shape_ptr(values.len(), values.as_ptr()) };
    let result = unsafe { PyArray1::borrow_from_array(&view, owner.clone().into_any()) };
    let _readonly = result.readwrite().make_nonwriteable();
    result
}

#[cfg(test)]
#[path = "interactions_tests.rs"]
mod tests;
