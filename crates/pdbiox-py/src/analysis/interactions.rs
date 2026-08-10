//! Spatial contacts and secondary structure in coarse Rust calls.

use crate::graph::PySpatialBackend;
use crate::structure::PyStructure;
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
    ) -> PyResult<Vec<PyContact>> {
        let structure = self.structure().clone();
        py.detach(move || pdbiox::analysis::atom_contacts(&structure, cutoff, backend.into()))
            .map(|contacts| contacts.into_iter().map(PyContact::from).collect())
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

impl From<pdbiox::analysis::Contact> for PyContact {
    fn from(value: pdbiox::analysis::Contact) -> Self {
        Self {
            first: value.first.get(),
            second: value.second.get(),
            distance: value.distance,
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
