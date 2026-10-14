//! Typed `NumPy` array arguments with one checked extraction.
//!
//! A `Bound<PyArray2<f32>>` parameter is checked by `PyO3` before the function
//! runs, and when the dtype is wrong the only thing it can say is that an
//! `ndarray` is not an `ndarray`. [`Array`] and [`Readonly`] perform the same
//! checks themselves and raise `MolframeTypeError` stating the dtype and rank
//! that were expected and the ones received. `PyO3` attaches the parameter name
//! to the exception as a note, so the traceback names the offending argument.
//!
//! No conversion is attempted: a `float64` array is rejected rather than copied
//! into `float32`, because the copy is the size of the input, narrows the values
//! silently and would make the cost of a call depend on the caller's dtype.
//! The caller opts in with `array.astype(numpy.float32)`.
//!
//! Every array parameter of the extension takes one of these two types, so the
//! behaviour and the wording are defined here and nowhere else.

use numpy::{
    Element, PyArray, PyArrayDescrMethods, PyReadonlyArray, PyUntypedArray, PyUntypedArrayMethods,
    dtype,
};
use pyo3::prelude::*;
use pyo3::pybacked::PyBackedStr;
use std::ops::Deref;

use numpy::ndarray::{Dimension, Ix1, Ix2, Ix3};

/// A shape-and-dtype-checked array argument, borrowed from the caller.
///
/// It holds a borrowed reference rather than an owned one, so it is `Copy` and
/// costs no reference-count traffic to pass around.
pub(crate) struct Array<'a, 'py, T: Element, D: Dimension>(Borrowed<'a, 'py, PyArray<T, D>>);

impl<T: Element, D: Dimension> Clone for Array<'_, '_, T, D> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<T: Element, D: Dimension> Copy for Array<'_, '_, T, D> {}

/// A one-dimensional [`Array`].
pub(crate) type Array1<'a, 'py, T> = Array<'a, 'py, T, Ix1>;
/// A two-dimensional [`Array`].
pub(crate) type Array2<'a, 'py, T> = Array<'a, 'py, T, Ix2>;
/// A three-dimensional [`Array`].
pub(crate) type Array3<'a, 'py, T> = Array<'a, 'py, T, Ix3>;

/// A checked array argument that is already read-locked for reading.
pub(crate) struct Readonly<'py, T: Element, D: Dimension>(PyReadonlyArray<'py, T, D>);

/// A read-locked one-dimensional [`Readonly`].
pub(crate) type Readonly1<'py, T> = Readonly<'py, T, Ix1>;
/// A read-locked two-dimensional [`Readonly`].
pub(crate) type Readonly2<'py, T> = Readonly<'py, T, Ix2>;

impl<'py, T: Element, D: Dimension> Deref for Array<'_, 'py, T, D> {
    type Target = Bound<'py, PyArray<T, D>>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<'py, T: Element, D: Dimension> Deref for Readonly<'py, T, D> {
    type Target = PyReadonlyArray<'py, T, D>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<'a, 'py, T: Element, D: Dimension> FromPyObject<'a, 'py> for Array<'a, 'py, T, D> {
    type Error = PyErr;

    fn extract(object: Borrowed<'a, 'py, PyAny>) -> PyResult<Self> {
        check::<T, D>(&object)?;
        // SAFETY: `check` has just verified the dtype and rank this type names.
        Ok(Self(unsafe { object.cast_unchecked() }))
    }
}

impl<'a, 'py, T: Element, D: Dimension> FromPyObject<'a, 'py> for Readonly<'py, T, D> {
    type Error = PyErr;

    fn extract(object: Borrowed<'a, 'py, PyAny>) -> PyResult<Self> {
        use numpy::PyArrayMethods;
        check::<T, D>(&object)?;
        // SAFETY: `check` has just verified the dtype and rank this type names.
        let array = unsafe { object.cast_unchecked::<PyArray<T, D>>() };
        Ok(Self(array.readonly()))
    }
}

/// Checks the object is an `ndarray` of dtype `T` and rank `D`.
fn check<T: Element, D: Dimension>(object: &Bound<'_, PyAny>) -> PyResult<()> {
    let py = object.py();
    let Ok(untyped) = object.cast::<PyUntypedArray>() else {
        let received = object.get_type().name()?;
        return Err(crate::error::type_error(not_an_array(&received.to_cow()?)));
    };
    let wanted = dtype::<T>(py);
    let actual = untyped.dtype();
    if !actual.is_equiv_to(&wanted) {
        return Err(crate::error::type_error(dtype_mismatch(
            &dtype_name(&wanted)?,
            &dtype_name(&actual)?,
        )));
    }
    match D::NDIM {
        Some(rank) if untyped.ndim() != rank => Err(crate::error::type_error(rank_mismatch(
            rank,
            untyped.ndim(),
        ))),
        _ => Ok(()),
    }
}

fn dtype_name(descriptor: &Bound<'_, numpy::PyArrayDescr>) -> PyResult<String> {
    let name = descriptor.as_any().getattr("name")?;
    Ok(name.extract::<PyBackedStr>()?.to_string())
}

fn not_an_array(received: &str) -> String {
    format!("expected a numpy.ndarray, received {received}")
}

fn dtype_mismatch(expected: &str, received: &str) -> String {
    format!(
        "expected a {expected} array, received dtype {received}; convert it once with \
         array.astype(numpy.{expected})"
    )
}

fn rank_mismatch(expected: usize, received: usize) -> String {
    format!("expected a {expected}-dimensional array, received {received} dimensions")
}

#[cfg(test)]
#[path = "array_tests.rs"]
mod tests;
