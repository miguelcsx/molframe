//! Coordinates as a `DLPack` tensor.
//!
//! `DLPack` has no read-only flag, so a consumer may legally write to a tensor it
//! imports. Exporting the structure's own storage would let it corrupt an
//! immutable snapshot, so every export is an independent copy.

use crate::bindings::PyStructure;
use molframe::interop::{DLManagedTensor, DlpackTensor};
use pyo3::exceptions::PyBufferError;
use pyo3::prelude::*;
use pyo3::{Bound, ffi};
use std::ffi::{CStr, c_void};

/// The name of a capsule that no consumer has taken.
const UNCONSUMED: &CStr = c"dltensor";
/// `DLPack`'s device code for host memory.
const CPU: i32 = 1;

/// The coordinates of a structure, as something `from_dlpack` can import.
///
/// `numpy.from_dlpack`, `torch.from_dlpack` and `jax.dlpack.from_dlpack` take
/// this object directly. Each import copies the `(atoms, 3)` `float32`
/// coordinates, so the consumer owns its tensor outright.
#[derive(Clone, Debug)]
#[pyclass(
    name = "CoordinateTensor",
    frozen,
    skip_from_py_object,
    module = "molframe.interop"
)]
pub(crate) struct PyCoordinateTensor {
    structure: PyStructure,
}

/// Frees a tensor whose capsule was destroyed without a consumer taking it.
///
/// A consumer that imports the tensor renames the capsule `used_dltensor` and
/// owns the deleter call from then on, so the name decides who frees: freeing a
/// tensor already taken would be a double free.
unsafe extern "C" fn release_unconsumed(capsule: *mut ffi::PyObject) {
    // SAFETY: CPython calls this with the capsule being destroyed. An exception
    // raised while it is destroyed belongs to the caller, so it is set aside
    // and put back whatever happens here.
    let pending = unsafe { ffi::PyErr_GetRaisedException() };
    // SAFETY: `capsule` is a live capsule; `PyCapsule_IsValid` only compares
    // its name and does not raise.
    if unsafe { ffi::PyCapsule_IsValid(capsule, UNCONSUMED.as_ptr()) } == 1 {
        // SAFETY: the name matched, so the pointer is the one stored below.
        let managed = unsafe { ffi::PyCapsule_GetPointer(capsule, UNCONSUMED.as_ptr()) };
        // SAFETY: still named `dltensor`, so nothing consumed the tensor and
        // this is the only release.
        unsafe { molframe::interop::release(managed.cast::<DLManagedTensor>()) };
    }
    // SAFETY: restores exactly what was set aside above (possibly nothing).
    unsafe { ffi::PyErr_SetRaisedException(pending) };
}

/// Wraps a managed tensor in the `dltensor` capsule a consumer imports.
fn capsule(py: Python<'_>, tensor: DlpackTensor) -> PyResult<Bound<'_, PyAny>> {
    let managed = tensor.into_raw();
    // SAFETY: `managed` is the pointer `into_raw` handed over and nothing else
    // frees it; on success the capsule owns it and `release_unconsumed` frees it
    // unless a consumer takes it first.
    let raw = unsafe {
        ffi::PyCapsule_New(
            managed.cast::<c_void>(),
            UNCONSUMED.as_ptr(),
            Some(release_unconsumed),
        )
    };
    if raw.is_null() {
        // SAFETY: the capsule was not made, so ownership never left this
        // function and the tensor has not been consumed.
        unsafe { molframe::interop::release(managed) };
        return Err(PyErr::fetch(py));
    }
    // SAFETY: `PyCapsule_New` returned a new strong reference.
    Ok(unsafe { Bound::from_owned_ptr(py, raw) })
}

#[pymethods]
impl PyCoordinateTensor {
    /// `(atoms, 3)`.
    #[getter]
    fn shape(&self) -> (usize, usize) {
        (self.structure.inner.atom_count() as usize, 3)
    }

    /// `"float32"`.
    #[getter]
    #[allow(clippy::unused_self)]
    const fn dtype(&self) -> &'static str {
        "float32"
    }

    /// A `dltensor` capsule over a private copy of the coordinates.
    #[pyo3(signature = (*, stream=None, max_version=None, dl_device=None, copy=None))]
    fn __dlpack__<'py>(
        &self,
        py: Python<'py>,
        stream: Option<&Bound<'py, PyAny>>,
        max_version: Option<(u32, u32)>,
        dl_device: Option<(i32, i32)>,
        copy: Option<bool>,
    ) -> PyResult<Bound<'py, PyAny>> {
        // Host memory has no stream to synchronise, and every export is a
        // copy, so `stream` and `copy` need no action; a legacy capsule is
        // what any `max_version` may be answered with.
        let _ = (stream, max_version, copy);
        if dl_device.is_some_and(|device| device != (CPU, 0)) {
            return Err(PyBufferError::new_err(
                "coordinates live in host memory; only the CPU device is available",
            ));
        }
        let structure = self.structure.inner.clone();
        let tensor = py
            .detach(move || DlpackTensor::coordinates(structure.engine()))
            .map_err(crate::error::kernel)?;
        capsule(py, tensor)
    }

    /// `(1, 0)`: host memory.
    #[allow(clippy::unused_self)]
    const fn __dlpack_device__(&self) -> (i32, i32) {
        (CPU, 0)
    }

    fn __repr__(&self) -> String {
        format!(
            "CoordinateTensor(shape=({}, 3), dtype=float32)",
            self.structure.inner.atom_count()
        )
    }
}

/// The coordinates of `structure` as a `DLPack` tensor source.
#[pyfunction]
pub(crate) fn coordinates(structure: &PyStructure) -> PyCoordinateTensor {
    PyCoordinateTensor {
        structure: structure.clone(),
    }
}
