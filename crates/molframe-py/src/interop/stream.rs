//! The Arrow C stream capsule every tabular object hands to a consumer.

use pyo3::prelude::*;
use pyo3::types::PyCapsule;

/// Wraps a native Arrow stream as the `arrow_array_stream` capsule that
/// `__arrow_c_stream__` returns.
///
/// The stream owns its producer and converts one batch per pull, so a consumer
/// may keep it, or the arrays it pulled, after every Python handle is gone.
pub(crate) fn arrow_capsule<E: std::fmt::Display>(
    py: Python<'_>,
    stream: Result<molframe::interop::ArrowStream, E>,
) -> PyResult<Bound<'_, PyCapsule>> {
    let stream = stream.map_err(|error| crate::error::internal(error.to_string()))?;
    PyCapsule::new_with_value(py, stream.into_ffi(), c"arrow_array_stream")
}
