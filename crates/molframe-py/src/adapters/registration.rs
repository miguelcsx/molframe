//! Registration of neutral transfer boundaries.

use pyo3::prelude::*;

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    super::bindings::register(module)
}
