//! Explicit capsule protocols shared by sibling Rust extensions.

use pyo3::prelude::*;
use pyo3::types::{PyAnyMethods, PyCapsule, PyCapsuleMethods};

const STRUCTURE_CAPSULE: &std::ffi::CStr = c"molframe.Structure";

/// Extracts a `molframe::Structure` through the stable Python capsule protocol.
///
/// # Errors
///
/// Returns a Python exception when the object exposes no capsule method or the
/// capsule fails its name check.
pub fn structure_from_python(object: &Bound<'_, PyAny>) -> PyResult<molframe::Structure> {
    let capsule = object.call_method0("_pdviewx_structure_capsule")?;
    let capsule = capsule.cast::<PyCapsule>()?;
    let pointer = capsule.pointer_checked(Some(STRUCTURE_CAPSULE))?;

    // SAFETY: `_pdviewx_structure_capsule` is the sole producer of this
    // capsule, its name check proves the protocol identity, and the capsule
    // remains alive while this reference is cloned. `Structure` is repr(C).
    Ok(unsafe { pointer.cast::<molframe::Structure>().as_ref() }.clone())
}
