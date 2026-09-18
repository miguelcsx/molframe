//! Contiguous immutable `NumPy` views over shared Rust trajectory buffers.

use numpy::ndarray::{ArrayView1, ArrayView2, ArrayView3};
use numpy::{PyArray1, PyArray2, PyArray3, PyArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::sync::Arc;

#[pyclass(frozen)]
pub(super) struct Float32Owner {
    pub(super) values: Arc<[f32]>,
}

#[pyclass(frozen)]
pub(super) struct Float64Owner {
    pub(super) values: Arc<[f64]>,
}

#[pyclass(frozen)]
pub(super) struct Int64Owner {
    pub(super) values: Arc<[i64]>,
}

pub(super) fn vectors(
    py: Python<'_>,
    values: Arc<[f32]>,
    frames: usize,
    atoms: usize,
) -> PyResult<Bound<'_, PyArray3<f32>>> {
    require_elements(values.len(), &[frames, atoms, 3])?;
    let owner = Bound::new(py, Float32Owner { values })?;
    let pointer = owner.borrow().values.as_ptr();
    // SAFETY: `require_elements` proved that the dimensions cover the complete
    // contiguous Arc allocation, whose element alignment matches `f32`.
    let view = unsafe { ArrayView3::from_shape_ptr((frames, atoms, 3), pointer) };
    // SAFETY: `owner` becomes the NumPy base, retaining that immutable Arc for
    // at least as long as the view; the result is made non-writeable below.
    let result = unsafe { PyArray3::borrow_from_array(&view, owner.into_any()) };
    let _readonly = result.readwrite().make_nonwriteable();
    Ok(result)
}

pub(super) fn scalars(py: Python<'_>, values: Arc<[f64]>) -> PyResult<Bound<'_, PyArray1<f64>>> {
    let owner = Bound::new(py, Float64Owner { values })?;
    let length = owner.borrow().values.len();
    let pointer = owner.borrow().values.as_ptr();
    // SAFETY: `length` is the exact element count of this aligned contiguous
    // `f64` Arc allocation.
    let view = unsafe { ArrayView1::from_shape_ptr(length, pointer) };
    // SAFETY: `owner` becomes the NumPy base and retains the immutable Arc;
    // the result is made non-writeable before it escapes.
    let result = unsafe { PyArray1::borrow_from_array(&view, owner.into_any()) };
    let _readonly = result.readwrite().make_nonwriteable();
    Ok(result)
}

pub(super) fn integers(py: Python<'_>, values: Arc<[i64]>) -> PyResult<Bound<'_, PyArray1<i64>>> {
    let owner = Bound::new(py, Int64Owner { values })?;
    let length = owner.borrow().values.len();
    let pointer = owner.borrow().values.as_ptr();
    // SAFETY: `length` is the exact element count of this aligned contiguous
    // `i64` Arc allocation.
    let view = unsafe { ArrayView1::from_shape_ptr(length, pointer) };
    // SAFETY: `owner` becomes the NumPy base and retains the immutable Arc;
    // the result is made non-writeable before it escapes.
    let result = unsafe { PyArray1::borrow_from_array(&view, owner.into_any()) };
    let _readonly = result.readwrite().make_nonwriteable();
    Ok(result)
}

pub(super) fn cells(
    py: Python<'_>,
    values: Arc<[f64]>,
    frames: usize,
) -> PyResult<Bound<'_, PyArray2<f64>>> {
    require_elements(values.len(), &[frames, 6])?;
    let owner = Bound::new(py, Float64Owner { values })?;
    let pointer = owner.borrow().values.as_ptr();
    // SAFETY: `require_elements` proved that the dimensions cover the complete
    // contiguous Arc allocation, whose element alignment matches `f64`.
    let view = unsafe { ArrayView2::from_shape_ptr((frames, 6), pointer) };
    // SAFETY: `owner` becomes the NumPy base, retaining that immutable Arc for
    // at least as long as the view; the result is made non-writeable below.
    let result = unsafe { PyArray2::borrow_from_array(&view, owner.into_any()) };
    let _readonly = result.readwrite().make_nonwriteable();
    Ok(result)
}

fn require_elements(actual: usize, dimensions: &[usize]) -> PyResult<()> {
    let expected = dimensions.iter().try_fold(1_usize, |product, dimension| {
        product.checked_mul(*dimension)
    });
    if expected == Some(actual) {
        Ok(())
    } else {
        Err(PyValueError::new_err(
            "trajectory buffer length does not match its declared shape",
        ))
    }
}

pub(super) fn frame_vectors(
    py: Python<'_>,
    lease: Arc<molframe::core::BatchLease<molframe::traj::TrajectoryBatch>>,
    column: u8,
) -> PyResult<Option<Bound<'_, PyArray2<f32>>>> {
    let owner = Bound::new(py, super::stream_frame::PyStreamFrame { lease })?;
    let (length, pointer) = {
        let data = owner.borrow();
        let frame = data
            .lease
            .batch()
            .timestep()
            .ok_or_else(|| PyValueError::new_err("frame is absent"))?;
        let values = match column {
            0 => Some(&frame.positions),
            1 => frame.velocities.as_ref(),
            _ => frame.forces.as_ref(),
        };
        let Some(values) = values else {
            return Ok(None);
        };
        (values.len(), values.as_ptr().cast::<f32>())
    };
    // SAFETY: the lease owns `length` contiguous triples. It cannot be recycled
    // while NumPy's base retains an Arc to the lease; mutation is forbidden.
    let view = unsafe { ArrayView2::from_shape_ptr((length, 3), pointer) };
    // SAFETY: owner retains the frame allocation and its memory reservation for
    // the entire array lifetime, including sliced arrays and shared references.
    let result = unsafe { PyArray2::borrow_from_array(&view, owner.into_any()) };
    let _readonly = result.readwrite().make_nonwriteable();
    Ok(Some(result))
}
