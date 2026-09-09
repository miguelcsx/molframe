//! Owned Python projections for non-crystallographic symmetry.

use super::types::{PyAffineTransform, PySymmetrySet};
use crate::cif_document::PyCifDocument;
use crate::structure::PyStructure;
use pyo3::prelude::*;

#[pyclass(name = "NcsCode", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyNcsCode {
    Given,
    Generate,
}

#[pyclass(name = "NcsOperator", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyNcsOperator {
    #[pyo3(get)]
    id: String,
    #[pyo3(get)]
    code: PyNcsCode,
    #[pyo3(get)]
    details: Option<String>,
    #[pyo3(get)]
    transform: PyAffineTransform,
}

#[pyclass(name = "NcsSet", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyNcsSet(pub(crate) pdbiox::xtal::NcsSet);

#[pymethods]
impl PyNcsSet {
    fn __len__(&self) -> usize {
        self.0.len()
    }
    #[getter]
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
    fn get(&self, id: &str) -> Option<PyNcsOperator> {
        self.0.get(id).cloned().map(Into::into)
    }
    fn operators(&self) -> Vec<PyNcsOperator> {
        self.0.operators().cloned().map(Into::into).collect()
    }
    fn generators(&self) -> Vec<PyNcsOperator> {
        self.0.generators().cloned().map(Into::into).collect()
    }
}

#[pyclass(name = "NcsView", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyNcsView {
    structure: pdbiox::Structure,
    set: pdbiox::xtal::NcsSet,
}

type NcsPosition = (u32, PyNcsOperator, Option<[f32; 3]>);

#[pymethods]
impl PyNcsView {
    #[getter]
    fn copy_count(&self) -> usize {
        self.set.generators().count()
    }

    fn atoms(&self) -> Vec<(u32, PyNcsOperator)> {
        self.set
            .generators()
            .flat_map(|operator| {
                (0..self.structure.atom_count()).map(move |atom| (atom, operator.clone().into()))
            })
            .collect()
    }

    #[pyo3(signature = (model=0))]
    fn positions(&self, model: usize) -> PyResult<Vec<NcsPosition>> {
        let model = model_index(model)?;
        let Some(positions) = self.structure.model_positions(model) else {
            return Ok(Vec::new());
        };
        Ok(self
            .set
            .generators()
            .flat_map(|operator| {
                (0..self.structure.atom_count()).map(move |atom| {
                    let position = positions
                        .get(atom as usize)
                        .copied()
                        .map(|position| operator.transform.apply(position));
                    (atom, operator.clone().into(), position)
                })
            })
            .collect())
    }
}

#[pyclass(name = "CrystalNeighbor", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyCrystalNeighbor {
    #[pyo3(get)]
    pub(crate) source_atom: u32,
    #[pyo3(get)]
    pub(crate) image_atom: u32,
    #[pyo3(get)]
    pub(crate) operation: usize,
    #[pyo3(get)]
    pub(crate) lattice: [i32; 3],
    #[pyo3(get)]
    pub(crate) distance_squared: f64,
}

#[pymethods]
impl PyStructure {
    fn symmetry_set(&self) -> Option<PySymmetrySet> {
        pdbiox::xtal::SymmetryExt::symmetry_set(self.structure())
            .cloned()
            .map(PySymmetrySet)
    }

    fn ncs_set(&self) -> Option<PyNcsSet> {
        pdbiox::xtal::NcsExt::ncs_set(self.structure())
            .cloned()
            .map(PyNcsSet)
    }

    fn ncs_generated(&self) -> Option<PyNcsView> {
        pdbiox::xtal::NcsExt::ncs_set(self.structure())
            .cloned()
            .map(|set| PyNcsView {
                structure: self.structure().clone(),
                set,
            })
    }
}

#[pyfunction]
pub(crate) fn lower_ncs(py: Python<'_>, document: &PyCifDocument) -> PyResult<PyNcsSet> {
    // Only the diagnostic conversion needs the interpreter.
    py.detach(|| pdbiox::xtal::lower_ncs(&document.inner))
        .map(PyNcsSet)
        .map_err(|findings| crate::errors::read_error(py, &findings))
}

fn model_index(model: usize) -> PyResult<pdbiox::ModelIndex> {
    u32::try_from(model)
        .map(pdbiox::ModelIndex::new)
        .map_err(|_| pyo3::exceptions::PyValueError::new_err("model index exceeds u32"))
}

impl From<pdbiox::xtal::NcsCode> for PyNcsCode {
    fn from(value: pdbiox::xtal::NcsCode) -> Self {
        match value {
            pdbiox::xtal::NcsCode::Given => Self::Given,
            pdbiox::xtal::NcsCode::Generate => Self::Generate,
        }
    }
}

impl From<pdbiox::xtal::NcsOperator> for PyNcsOperator {
    fn from(value: pdbiox::xtal::NcsOperator) -> Self {
        Self {
            id: value.id.to_string(),
            code: value.code.into(),
            details: value.details.map(|value| value.to_string()),
            transform: PyAffineTransform(value.transform),
        }
    }
}

impl From<pdbiox::xtal::CrystalNeighbor> for PyCrystalNeighbor {
    fn from(value: pdbiox::xtal::CrystalNeighbor) -> Self {
        Self {
            source_atom: value.source_atom.get(),
            image_atom: value.image_atom.get(),
            operation: value.operation,
            lattice: value.lattice,
            distance_squared: value.distance_squared,
        }
    }
}
