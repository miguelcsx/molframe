//! Borrowing native structures out of Python objects.

use super::bindings;
use bindings::PyStructure;
use pyo3::prelude::*;

/// Borrows the native snapshot retained by a Python `molframe.Structure`.
///
/// Cloning the returned value only increments its shared storage ownership;
/// coordinate columns are not copied.
///
/// # Errors
///
/// Returns Python's type error when `object` is not a native `MolFrame`
/// structure.
pub fn structure_from_python(object: &Bound<'_, PyAny>) -> PyResult<molframe::Structure> {
    Ok(object
        .extract::<PyRef<'_, PyStructure>>()
        .map(|structure| structure.inner.clone())?)
}
