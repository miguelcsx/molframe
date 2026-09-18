//! File-backed declarative functional-evaluation specifications.

use super::errors::specification_error;
use super::specification::PyMotif;
use super::verdict::PyVerdictProfile;
use pyo3::prelude::*;
use std::path::PathBuf;

#[pyclass(name = "EvaluationSpecification", frozen, skip_from_py_object)]
#[derive(Debug)]
pub(crate) struct PyEvaluationSpecification {
    inner: molframe::fx::EvaluationSpecification,
}

#[pymethods]
impl PyEvaluationSpecification {
    #[getter]
    fn motif(&self) -> PyMotif {
        PyMotif(self.inner.motif.clone())
    }

    #[getter]
    fn profile(&self) -> PyVerdictProfile {
        PyVerdictProfile(self.inner.profile.clone())
    }
}

#[pyfunction]
pub(crate) fn read_evaluation_specification(
    py: Python<'_>,
    path: PathBuf,
) -> PyResult<PyEvaluationSpecification> {
    py.detach(move || molframe::fx::read_evaluation_specification(path))
        .map(|inner| PyEvaluationSpecification { inner })
        .map_err(specification_error)
}
