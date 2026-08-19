//! Native graph-equivalence enumeration at the Python boundary.

use super::components::PyComponent;
use ::pdbiox;
use pyo3::create_exception;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

create_exception!(_native, AutomorphismLimit, PyValueError);

#[pyfunction]
pub(crate) fn equivalence_classes(
    component: &PyComponent,
) -> super::components::PyEquivalenceClasses {
    super::components::PyEquivalenceClasses(pdbiox::equivalence_classes(&component.0))
}

#[pyfunction]
pub(crate) fn automorphisms(component: &PyComponent, limit: usize) -> PyResult<Vec<Vec<u32>>> {
    pdbiox::automorphisms(&component.0, limit)
        .map(|mappings| mappings.into_iter().map(<[u32]>::into_vec).collect())
        .map_err(|error| AutomorphismLimit::new_err(error.to_string()))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add(
        "AutomorphismLimit",
        module.py().get_type::<AutomorphismLimit>(),
    )?;
    module.add_function(wrap_pyfunction!(equivalence_classes, module)?)?;
    module.add_function(wrap_pyfunction!(automorphisms, module)?)?;
    Ok(())
}
