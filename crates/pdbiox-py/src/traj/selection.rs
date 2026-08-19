//! Fixed-topology selection plans for trajectory frames.

use super::reader_types::PyTimestep;
use crate::graph::PySpatialBackend;
use crate::query::{PyAnalysisPolicy, PyEvaluation, PyQuery, PySelection, groups_from_python};
use crate::structure::PyStructure;
use pyo3::prelude::*;
use std::collections::BTreeMap;
use std::sync::Arc;

/// A compiled query that reuses topology planning across frame evaluations.
#[pyclass(name = "UpdatingSelection", frozen, skip_from_py_object)]
pub(crate) struct PyUpdatingSelection {
    inner: Arc<pdbiox::traj::UpdatingSelection>,
}

#[pymethods]
impl PyUpdatingSelection {
    #[new]
    #[pyo3(signature = (topology, query, *, policy=None, groups=None, backend=PySpatialBackend::Auto))]
    fn new(
        topology: &PyStructure,
        query: &PyQuery,
        policy: Option<&PyAnalysisPolicy>,
        groups: Option<BTreeMap<String, PySelection>>,
        backend: PySpatialBackend,
    ) -> Self {
        let policy =
            policy.map_or_else(pdbiox::AnalysisPolicy::default, |value| value.inner.clone());
        let groups = groups_from_python(groups);
        let inner = pdbiox::traj::UpdatingSelection::new(
            topology.structure().clone(),
            &query.inner,
            policy,
            groups,
            backend.into(),
        );
        Self {
            inner: Arc::new(inner),
        }
    }

    fn evaluate(&self, py: Python<'_>, timestep: &PyTimestep) -> PyResult<PyEvaluation> {
        let frame = timestep.clone().try_into()?;
        let selection = Arc::clone(&self.inner);
        py.detach(move || selection.evaluate(&frame))
            .map(PyEvaluation::from)
            .map_err(|error| crate::errors::UpdatingSelectionError::new_err(error.to_string()))
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyUpdatingSelection>()?;
    Ok(())
}
