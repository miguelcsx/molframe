//! Explicit coordinate-frame decisions for functional measurements.

use super::mapping::PyMappedMotif;
use crate::geometry::PyRigid;
use pyo3::prelude::*;

#[pyclass(name = "AlignmentKind", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyAlignmentKind {
    NotRequired,
    CallerSupplied,
}

impl From<PyAlignmentKind> for molframe::fx::AlignmentKind {
    fn from(value: PyAlignmentKind) -> Self {
        match value {
            PyAlignmentKind::NotRequired => Self::NotRequired,
            PyAlignmentKind::CallerSupplied => Self::CallerSupplied,
        }
    }
}

impl From<molframe::fx::AlignmentKind> for PyAlignmentKind {
    fn from(value: molframe::fx::AlignmentKind) -> Self {
        match value {
            molframe::fx::AlignmentKind::NotRequired => Self::NotRequired,
            molframe::fx::AlignmentKind::CallerSupplied => Self::CallerSupplied,
        }
    }
}

#[pyclass(name = "AlignedMotif", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAlignedMotif(pub(crate) molframe::fx::AlignedMotif);

#[pymethods]
impl PyAlignedMotif {
    #[new]
    fn new(mapping: &PyMappedMotif, transform: &PyRigid, kind: PyAlignmentKind) -> Self {
        Self(molframe::fx::AlignedMotif {
            mapping: mapping.0.clone(),
            transform: transform.0,
            kind: kind.into(),
        })
    }

    #[getter]
    fn mapping(&self) -> PyMappedMotif {
        PyMappedMotif(self.0.mapping.clone())
    }

    #[getter]
    fn transform(&self) -> PyRigid {
        PyRigid(self.0.transform)
    }

    #[getter]
    fn kind(&self) -> PyAlignmentKind {
        self.0.kind.into()
    }
}

#[pyfunction]
pub(crate) fn align_intrinsic(py: Python<'_>, mapping: &PyMappedMotif) -> PyAlignedMotif {
    py.detach(move || -> PyAlignedMotif {
        PyAlignedMotif(molframe::fx::align_intrinsic(mapping.0.clone()))
    })
}

#[pyfunction]
pub(crate) fn align_with_transform(
    py: Python<'_>,
    mapping: &PyMappedMotif,
    transform: &PyRigid,
) -> PyAlignedMotif {
    py.detach(move || -> PyAlignedMotif {
        PyAlignedMotif(molframe::fx::align_with_transform(
            mapping.0.clone(),
            transform.0,
        ))
    })
}
