//! Versioned native structure sharing between independently loaded extensions.
//!
//! The capsule contains only C-layout records, primitive slices and one C-ABI
//! query callback. Rust-owned values never cross the extension boundary. The
//! consumer retains the capsule, which in turn retains the immutable structure
//! snapshot and every exported buffer.

use pyo3::prelude::*;
use pyo3::types::{PyAny, PyCapsule, PyCapsuleMethods};
use std::ffi::{CStr, c_char};
use std::fmt;
use std::slice;
use std::sync::OnceLock;

#[cfg(test)]
#[path = "native_source_tests.rs"]
mod tests;

const CAPSULE_NAME: &CStr = c"molframe.StructureSource.v1";
const ABI_VERSION: u32 = 1;

/// Compact renderer-facing atom metadata.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[repr(C)]
pub struct NativeAtom {
    /// Atomic number, or zero when unknown.
    pub element: u16,
    /// Owning residue row.
    pub residue: u32,
}

/// Compact covalent-bond endpoints.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[repr(C)]
pub struct NativeBond {
    /// First atom row.
    pub first: u32,
    /// Second atom row.
    pub second: u32,
}

/// Owned compact topology copied once from the parser's compressed columns.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NativeTopology {
    /// Atom metadata in coordinate order.
    pub atoms: Vec<NativeAtom>,
    /// Residue-to-atom start offsets including the final end offset.
    pub residue_atom_start: Vec<u32>,
    /// Chain-to-residue start offsets including the final end offset.
    pub chain_residue_start: Vec<u32>,
    /// Model-to-chain start offsets including the final end offset.
    pub model_chain_start: Vec<u32>,
    /// Covalent bonds.
    pub bonds: Vec<NativeBond>,
}

type SelectFn =
    unsafe extern "C" fn(*const NativeSourceV1, *const c_char, usize, *mut u32, usize) -> i64;
type EncodeFn = unsafe extern "C" fn(*const NativeSourceV1, *mut u8, usize) -> i64;

#[repr(C)]
struct NativeSourceV1 {
    abi_version: u32,
    coordinate_generation: u64,
    coordinates: *const f32,
    coordinate_rows: usize,
    atoms: *const NativeAtom,
    atom_count: usize,
    residue_atom_start: *const u32,
    residue_atom_start_len: usize,
    chain_residue_start: *const u32,
    chain_residue_start_len: usize,
    model_chain_start: *const u32,
    model_chain_start_len: usize,
    bonds: *const NativeBond,
    bond_count: usize,
    select: SelectFn,
    encode_bcif: EncodeFn,
}

// The raw pointers refer only to immutable buffers retained by the same
// capsule. Python owns capsule destruction, and access requires an attached
// interpreter through `NativeStructureSource`.
unsafe impl Send for NativeSourceV1 {}

#[repr(C)]
struct NativeCapsule {
    api: NativeSourceV1,
    structure: molframe::Structure,
    topology: NativeTopology,
    encoded_bcif: OnceLock<Result<Vec<u8>, ()>>,
}

unsafe impl Send for NativeCapsule {}

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
        let capsule = object.call_method0("_molframe_source_v1")?;
        let capsule = capsule.cast_into::<PyCapsule>()?;
        let pointer = capsule.pointer_checked(Some(CAPSULE_NAME))?;
        let api = pointer.as_ptr() as usize;
        // SAFETY: the checked capsule name identifies a `NativeSourceV1`, its
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
    /// # Errors
    ///
    /// Returns a value error for malformed or unsupported queries.
    pub fn select(&self, source: &str) -> PyResult<Vec<u32>> {
        let capacity = self.coordinates().len();
        let mut rows = vec![0u32; capacity];
        // SAFETY: the callback belongs to the retained producer, the query and
        // output buffers remain valid for the call, and capacity is exact.
        let count = unsafe {
            let api = api_ref(self.api);
            (api.select)(
                api,
                source.as_ptr().cast::<c_char>(),
                source.len(),
                rows.as_mut_ptr(),
                rows.len(),
            )
        };
        let count = usize::try_from(count).map_err(|_| {
            pyo3::exceptions::PyValueError::new_err("MolFrame query evaluation failed")
        })?;
        rows.truncate(count);
        Ok(rows)
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

unsafe fn api_ref<'a>(address: usize) -> &'a NativeSourceV1 {
    // SAFETY: callers obtain this address only from a checked, retained capsule.
    unsafe { &*(address as *const NativeSourceV1) }
}

pub(crate) fn capsule<'py>(
    py: Python<'py>,
    structure: &molframe::Structure,
) -> PyResult<Bound<'py, PyCapsule>> {
    let topology = compact_topology(structure);
    let coordinates = structure.coordinates();
    let api = NativeSourceV1 {
        abi_version: ABI_VERSION,
        coordinate_generation: structure.engine().generation().get(),
        coordinates: coordinates.as_ptr().cast::<f32>(),
        coordinate_rows: coordinates.len(),
        atoms: topology.atoms.as_ptr(),
        atom_count: topology.atoms.len(),
        residue_atom_start: topology.residue_atom_start.as_ptr(),
        residue_atom_start_len: topology.residue_atom_start.len(),
        chain_residue_start: topology.chain_residue_start.as_ptr(),
        chain_residue_start_len: topology.chain_residue_start.len(),
        model_chain_start: topology.model_chain_start.as_ptr(),
        model_chain_start_len: topology.model_chain_start.len(),
        bonds: topology.bonds.as_ptr(),
        bond_count: topology.bonds.len(),
        select,
        encode_bcif,
    };
    PyCapsule::new_with_value(
        py,
        NativeCapsule {
            api,
            structure: structure.clone(),
            topology,
            encoded_bcif: OnceLock::new(),
        },
        CAPSULE_NAME,
    )
}

unsafe extern "C" fn encode_bcif(
    source: *const NativeSourceV1,
    output: *mut u8,
    capacity: usize,
) -> i64 {
    // SAFETY: the API is the first field of its retained capsule allocation.
    let capsule = unsafe { &*source.cast::<NativeCapsule>() };
    let encoded = capsule
        .encoded_bcif
        .get_or_init(|| molframe::write_bcif(&capsule.structure).map_err(|_| ()));
    let Ok(encoded) = encoded else {
        return -1;
    };
    let Ok(length) = i64::try_from(encoded.len()) else {
        return -1;
    };
    if output.is_null() || capacity == 0 {
        return length;
    }
    if capacity < encoded.len() {
        return -1;
    }
    // SAFETY: the caller supplies at least `encoded.len()` writable bytes.
    unsafe { std::ptr::copy_nonoverlapping(encoded.as_ptr(), output, encoded.len()) };
    length
}

fn compact_topology(structure: &molframe::Structure) -> NativeTopology {
    let data = structure.engine().data();
    let atoms = data
        .atoms()
        .map(|atom| NativeAtom {
            element: atom
                .element()
                .map_or(0, |element| u16::from(element.atomic_number())),
            residue: atom.residue().map_or(0, |residue| residue.index().get()),
        })
        .collect();
    let residue_atom_start = offsets(
        data.residues()
            .map(|residue| residue.atoms().map(|atom| atom.index().get())),
    );
    let chain_residue_start = offsets(
        data.chains()
            .map(|chain| chain.residues().map(|residue| residue.index().get())),
    );
    let model_chain_start = offsets(
        data.models()
            .map(|model| model.chains().map(|chain| chain.index().get())),
    );
    let bonds = structure
        .bonds()
        .iter()
        .map(|bond| NativeBond {
            first: bond.atom_a.get(),
            second: bond.atom_b.get(),
        })
        .collect();
    NativeTopology {
        atoms,
        residue_atom_start,
        chain_residue_start,
        model_chain_start,
        bonds,
    }
}

fn offsets<I, R>(rows: I) -> Vec<u32>
where
    I: Iterator<Item = R>,
    R: Iterator<Item = u32>,
{
    let mut starts = Vec::new();
    let mut end = 0;
    for row in rows {
        let mut row = row.peekable();
        let mut start = end;
        if let Some(value) = row.peek().copied() {
            start = value;
        }
        starts.push(start);
        end = row.last().map_or(start, |value| value.saturating_add(1));
    }
    starts.push(end);
    starts
}

unsafe extern "C" fn select(
    api: *const NativeSourceV1,
    source: *const c_char,
    source_len: usize,
    output: *mut u32,
    output_len: usize,
) -> i64 {
    if api.is_null() || source.is_null() || output.is_null() {
        return -1;
    }
    // SAFETY: the consumer passes back the API pointer supplied by this
    // producer; `NativeSourceV1` is the first field of `NativeCapsule`.
    let capsule = unsafe { &*api.cast::<NativeCapsule>() };
    // SAFETY: the caller promises a readable query buffer of `source_len`.
    let bytes = unsafe { slice::from_raw_parts(source.cast::<u8>(), source_len) };
    let Ok(source) = std::str::from_utf8(bytes) else {
        return -2;
    };
    let Ok(selection) = capsule
        .structure
        .select(source, &molframe::AnalysisPolicy::default())
    else {
        return -3;
    };
    let count = match usize::try_from(selection.len()) {
        Ok(value) if value <= output_len => value,
        _ => return -4,
    };
    // SAFETY: the consumer supplies a writable buffer of `output_len` rows.
    let output = unsafe { slice::from_raw_parts_mut(output, output_len) };
    for (slot, atom) in output.iter_mut().zip(selection.atoms()) {
        *slot = atom.index().get();
    }
    let mut encoded = -5;
    if let Ok(value) = i64::try_from(count) {
        encoded = value;
    }
    encoded
}
