//! Mutable-safe `DLPack` protocol export delegated to Rust.

use crate::structure::PyStructure;
use pyo3::exceptions::{PyBufferError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyCapsule;
use std::ffi::c_void;
use std::ptr::NonNull;

const CPU_DEVICE_TYPE: i32 = 1;
const CPU_DEVICE_ID: i32 = 0;

#[pymethods]
impl PyStructure {
    /// Returns the standard CPU device tuple used by array consumers.
    #[staticmethod]
    fn __dlpack_device__() -> (i32, i32) {
        (CPU_DEVICE_TYPE, CPU_DEVICE_ID)
    }

    /// Publishes a one-shot coordinate tensor capsule with independent storage.
    #[pyo3(signature = (stream=None, *, max_version=None, dl_device=None, copy=None))]
    fn __dlpack__<'py>(
        &self,
        py: Python<'py>,
        stream: Option<&Bound<'py, PyAny>>,
        max_version: Option<(u32, u32)>,
        dl_device: Option<(i32, i32)>,
        copy: Option<bool>,
    ) -> PyResult<Bound<'py, PyCapsule>> {
        let _ = (stream, max_version);
        if copy == Some(false) {
            return Err(PyBufferError::new_err(
                "copy=False is unavailable because DLPack cannot mark immutable storage",
            ));
        }
        if dl_device.is_some_and(|device| device != (CPU_DEVICE_TYPE, CPU_DEVICE_ID)) {
            return Err(PyValueError::new_err(
                "pdbiox structure coordinates are available on CPU device 0",
            ));
        }
        let tensor = pdbiox::DlpackTensor::coordinates(self.structure())
            .map_err(|error| PyBufferError::new_err(error.to_string()))?;
        let pointer = NonNull::new(tensor.into_raw().cast::<c_void>())
            .ok_or_else(|| PyBufferError::new_err("DLPack ownership transfer failed"))?;
        // SAFETY: the Rust producer transferred the managed tensor to this
        // capsule. Its DLPack deleter owns both metadata and independent tensor
        // storage, so a mutable consumer cannot alter a Structure snapshot.
        unsafe {
            PyCapsule::new_with_pointer_and_destructor(
                py,
                pointer,
                c"dltensor",
                Some(delete_unused_capsule),
            )
        }
    }
}

/// Releases a `DLPack` producer allocation when its capsule was never consumed.
///
/// # Safety
///
/// `CPython` must call this only with the live capsule that owns the producer
/// pointer; all raw access below is additionally gated by its protocol name.
unsafe extern "C" fn delete_unused_capsule(capsule: *mut pyo3::ffi::PyObject) {
    // A consumer renames an accepted capsule to `used_dltensor`, so only an
    // unconsumed capsule still validates under the producer name.
    // SAFETY: the destructor contract provides a live capsule pointer and the
    // protocol name is a static, null-terminated C string.
    let valid = unsafe { pyo3::ffi::PyCapsule_IsValid(capsule, c"dltensor".as_ptr()) } != 0;
    if !valid {
        return;
    }
    // SAFETY: `PyCapsule_IsValid` above proved the live capsule and exact name.
    let pointer = unsafe { pyo3::ffi::PyCapsule_GetPointer(capsule, c"dltensor".as_ptr()) }
        .cast::<pdbiox::DLManagedTensor>();
    let Some(tensor) = NonNull::new(pointer) else {
        return;
    };
    // SAFETY: validation above proves this is the producer-owned managed
    // tensor. DLPack requires its published deleter to release it exactly once.
    if let Some(deleter) = unsafe { tensor.as_ref() }.deleter {
        // SAFETY: this unconsumed capsule still owns the one-shot producer
        // pointer, so invoking its published deleter satisfies DLPack transfer.
        unsafe { deleter(tensor.as_ptr()) };
    }
}

#[cfg(test)]
#[path = "dlpack_tests.rs"]
mod tests;
