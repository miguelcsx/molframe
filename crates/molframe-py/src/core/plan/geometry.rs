//! Declarative geometry operations over borrowed `NumPy` buffers.

#[path = "geometry/classes.rs"]
mod classes;
#[path = "geometry/model.rs"]
mod model;
#[path = "geometry/native.rs"]
mod native;
#[path = "geometry/serialization.rs"]
mod serialization;
#[path = "geometry/values.rs"]
mod values;

pub(crate) use model::PyGeometryOperation;
pub(super) use native::{
    eigen_from_dict, position_from_dict, principal_from_dict, serialized_request, to_native,
    typed_request, validate_operation, weighted_from_dict,
};
pub(super) use serialization::{explain_operation, operation_to_dict};
pub(super) use values::{execute_operation, value_to_python};

use pyo3::prelude::*;

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    classes::register(module)
}
