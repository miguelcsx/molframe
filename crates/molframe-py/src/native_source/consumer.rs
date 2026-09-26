//! The safe view a consuming extension holds over a producer's capsule.

use super::abi::{ABI_VERSION, CAPSULE_NAME, NativeSourceV2, NativeTopology};
use super::producer::SELECT_NULL_ARGUMENT;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyCapsule, PyCapsuleMethods};
use std::ffi::c_char;
use std::fmt;
use std::slice;

/// Safe retained view of a `molframe.Structure` owned by another extension.
pub struct NativeStructureSource {
    capsule: Py<PyCapsule>,
    api: usize,
}

impl fmt::Debug for NativeStructureSource {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("NativeStructureSource")
            .field("coordinate_rows", &self.coordinates().len())
            .field("coordinate_generation", &self.coordinate_generation())
            .finish_non_exhaustive()
    }
}

impl Clone for NativeStructureSource {
    fn clone(&self) -> Self {
        Python::attach(|py| Self {
            capsule: self.capsule.clone_ref(py),
            api: self.api,
        })
    }
}

impl NativeStructureSource {
    /// Imports the versioned native provider exposed by a Python structure.
    ///
    /// # Errors
    ///
    /// Returns a Python exception when the object is not a compatible
    /// `molframe.Structure` or exposes a different ABI version.
    pub fn from_python(object: &Bound<'_, PyAny>) -> PyResult<Self> {
        let capsule = object.call_method0("_molframe_source_v2")?;
        let capsule = capsule.cast_into::<PyCapsule>()?;
        let pointer = capsule.pointer_checked(Some(CAPSULE_NAME))?;
        let api = pointer.as_ptr() as usize;
        // SAFETY: the checked capsule name identifies a `NativeSourceV2`, its
        // first field is the ABI version, and the producer retains its storage.
        let version = unsafe { api_ref(api).abi_version };
        if version != ABI_VERSION {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "unsupported MolFrame native source ABI {version}"
            )));
        }
        Ok(Self {
            capsule: capsule.unbind(),
            api,
        })
    }

    /// Parser-owned coordinate rows without an intermediate allocation.
    #[must_use]
    pub fn coordinates(&self) -> &[[f32; 3]] {
        // SAFETY: the versioned capsule retains a contiguous coordinate block
        // of `coordinate_rows * 3` finite-or-NaN `f32` values for this handle.
        unsafe {
            let api = api_ref(self.api);
            slice::from_raw_parts(api.coordinates.cast::<[f32; 3]>(), api.coordinate_rows)
        }
    }

    /// Coordinate generation used to invalidate spatial caches.
    #[must_use]
    pub fn coordinate_generation(&self) -> u64 {
        // SAFETY: `self.api` remains valid while `self.capsule` is retained.
        unsafe { api_ref(self.api).coordinate_generation }
    }

    /// Copies compact topology columns once for renderer-side dense access.
    #[must_use]
    pub fn topology(&self) -> NativeTopology {
        // SAFETY: every pointer and length belongs to an immutable vector
        // retained by the checked capsule for the lifetime of this handle.
        unsafe {
            let api = api_ref(self.api);
            NativeTopology {
                atoms: slice::from_raw_parts(api.atoms, api.atom_count).to_vec(),
                residue_atom_start: slice::from_raw_parts(
                    api.residue_atom_start,
                    api.residue_atom_start_len,
                )
                .to_vec(),
                chain_residue_start: slice::from_raw_parts(
                    api.chain_residue_start,
                    api.chain_residue_start_len,
                )
                .to_vec(),
                model_chain_start: slice::from_raw_parts(
                    api.model_chain_start,
                    api.model_chain_start_len,
                )
                .to_vec(),
                bonds: slice::from_raw_parts(api.bonds, api.bond_count).to_vec(),
            }
        }
    }

    /// Evaluates `MolFrame`'s canonical query and returns atom rows.
    ///
    /// The producer is first asked how many rows the query selects, and then
    /// to copy exactly that many, so a handful of selected atoms costs a
    /// handful of rows rather than one row per atom in the structure. The
    /// producer keeps the evaluated rows between the two calls, so the query is
    /// evaluated once. A producer that predates the size request refuses it;
    /// the rows are then copied into a buffer sized for the whole structure,
    /// which is always large enough, and trimmed afterwards.
    ///
    /// # Errors
    ///
    /// Returns a value error for malformed or unsupported queries.
    pub fn select(&self, source: &str) -> PyResult<Vec<u32>> {
        let failed = || pyo3::exceptions::PyValueError::new_err("MolFrame query evaluation failed");
        // SAFETY: the callback belongs to the retained producer and the query
        // buffer is valid for the call. A null output asks only for the count.
        let needed = unsafe { self.select_into(source, std::ptr::null_mut(), 0) };
        let capacity = match usize::try_from(needed) {
            Ok(count) => count,
            Err(_) if needed == SELECT_NULL_ARGUMENT => self.coordinates().len(),
            Err(_) => return Err(failed()),
        };
        let mut rows = vec![0u32; capacity];
        // SAFETY: `rows` is writable for `capacity` rows for the whole call.
        let count = unsafe { self.select_into(source, rows.as_mut_ptr(), rows.len()) };
        let count = usize::try_from(count).map_err(|_| failed())?;
        rows.truncate(count);
        rows.shrink_to_fit();
        Ok(rows)
    }

    /// Calls the producer's selection callback once.
    ///
    /// # Safety
    ///
    /// `output` must be null with `capacity` zero, or writable for `capacity`
    /// rows for the duration of the call.
    unsafe fn select_into(&self, source: &str, output: *mut u32, capacity: usize) -> i64 {
        // SAFETY: `self.api` remains valid while `self.capsule` is retained,
        // and the caller upholds the output contract.
        unsafe {
            let api = api_ref(self.api);
            (api.select)(
                api,
                source.as_ptr().cast::<c_char>(),
                source.len(),
                output,
                capacity,
            )
        }
    }

    /// Lazily encodes and caches browser-ready BCIF bytes in the producer.
    ///
    /// # Errors
    ///
    /// Returns a value error when the source cannot be encoded.
    pub fn encode_bcif(&self) -> PyResult<Vec<u8>> {
        // SAFETY: the retained checked capsule owns the callback and source.
        let size = unsafe {
            let api = api_ref(self.api);
            (api.encode_bcif)(api, std::ptr::null_mut(), 0)
        };
        let size = usize::try_from(size).map_err(|_| {
            pyo3::exceptions::PyValueError::new_err("MolFrame BCIF encoding failed")
        })?;
        let mut bytes = vec![0; size];
        // SAFETY: `bytes` has the exact capacity returned by the first call.
        let written = unsafe {
            let api = api_ref(self.api);
            (api.encode_bcif)(api, bytes.as_mut_ptr(), bytes.len())
        };
        if usize::try_from(written).ok() != Some(size) {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "MolFrame BCIF encoding failed",
            ));
        }
        Ok(bytes)
    }
}

unsafe fn api_ref<'a>(address: usize) -> &'a NativeSourceV2 {
    // SAFETY: callers obtain this address only from a checked, retained capsule.
    unsafe { &*(address as *const NativeSourceV2) }
}
