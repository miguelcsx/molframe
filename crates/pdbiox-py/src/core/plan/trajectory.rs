//! Declarative trajectory analyses over borrowed `NumPy` frame buffers.

#[path = "trajectory/classes.rs"]
mod classes;
#[path = "trajectory/model.rs"]
mod model;
#[path = "trajectory/native.rs"]
mod native;
#[path = "trajectory/serialization.rs"]
mod serialization;
#[path = "trajectory/values.rs"]
mod values;

pub(crate) use model::PyTrajectoryOperation;
pub(super) use native::{serialized_request, to_native, typed_request};
pub(super) use serialization::{explain_operation, operation_to_dict};
pub(super) use values::{execute_operation, value_to_python};

use pyo3::prelude::*;

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    classes::register(module)?;
    values::register(module)
}
