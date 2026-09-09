//! Immutable `NumPy` ownership for execution-charged f64 results.

use numpy::{PyArray1, PyArrayMethods, ndarray::ArrayView1};
use pyo3::prelude::*;

#[pyclass(frozen)]
struct Retained64Owner {
    values: pdbiox::core::execution::Retained<Vec<f64>>,
}

pub(crate) fn retained_scalars(
    py: Python<'_>,
    values: pdbiox::core::execution::Retained<Vec<f64>>,
) -> PyResult<Bound<'_, PyArray1<f64>>> {
    let owner = Bound::new(py, Retained64Owner { values })?;
    let (length, pointer) = {
        let data = owner.borrow();
        (data.values.len(), data.values.as_ptr())
    };
    // SAFETY: the immutable Vec contains exactly `length` contiguous f64 values.
    let view = unsafe { ArrayView1::from_shape_ptr(length, pointer) };
    // SAFETY: NumPy retains `owner`, including the Vec and its reservation, as
    // its base. Neither can be mutated or dropped while any array view lives.
    let result = unsafe { PyArray1::borrow_from_array(&view, owner.into_any()) };
    let _readonly = result.readwrite().make_nonwriteable();
    Ok(result)
}
