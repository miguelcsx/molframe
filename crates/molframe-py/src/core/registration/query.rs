//! Registration for compiled query and typed builder objects.

use crate::query::{
    PyAltlocPolicy, PyAnalysisPolicy, PyAssemblyChoice, PyEvaluation, PyLogicalPlan,
    PyMissingPolicy, PyModelChoice, PyNamespace, PyPhysicalQuery, PyQuery, PySelection,
};
use pyo3::prelude::*;

pub(super) fn register_query(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyNamespace>()?;
    module.add_class::<PyMissingPolicy>()?;
    module.add_class::<PyModelChoice>()?;
    module.add_class::<PyAltlocPolicy>()?;
    module.add_class::<PyAssemblyChoice>()?;
    module.add_class::<PyAnalysisPolicy>()?;
    module.add_class::<PySelection>()?;
    module.add_class::<PyEvaluation>()?;
    module.add_class::<PyQuery>()?;
    module.add_class::<PyLogicalPlan>()?;
    module.add_class::<PyPhysicalQuery>()?;
    crate::query::register(module)
}
