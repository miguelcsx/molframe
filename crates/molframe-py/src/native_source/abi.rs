//! The C layout shared by a structure's producer and its consumers.

use std::ffi::{CStr, c_char};
use std::sync::{Mutex, OnceLock};

pub(super) const CAPSULE_NAME: &CStr = c"molframe.StructureSource.v2";
pub(super) const ABI_VERSION: u32 = 2;

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
    /// One when the source assigns aromatic order.
    pub aromatic: u8,
    pub(super) reserved: [u8; 3],
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

pub(super) type SelectFn =
    unsafe extern "C" fn(*const NativeSourceV2, *const c_char, usize, *mut u32, usize) -> i64;
pub(super) type EncodeFn = unsafe extern "C" fn(*const NativeSourceV2, *mut u8, usize) -> i64;

#[repr(C)]
pub(super) struct NativeSourceV2 {
    pub(super) abi_version: u32,
    pub(super) coordinate_generation: u64,
    pub(super) coordinates: *const f32,
    pub(super) coordinate_rows: usize,
    pub(super) atoms: *const NativeAtom,
    pub(super) atom_count: usize,
    pub(super) residue_atom_start: *const u32,
    pub(super) residue_atom_start_len: usize,
    pub(super) chain_residue_start: *const u32,
    pub(super) chain_residue_start_len: usize,
    pub(super) model_chain_start: *const u32,
    pub(super) model_chain_start_len: usize,
    pub(super) bonds: *const NativeBond,
    pub(super) bond_count: usize,
    pub(super) select: SelectFn,
    pub(super) encode_bcif: EncodeFn,
}

// The raw pointers refer only to immutable buffers retained by the same
// capsule. Python owns capsule destruction, and access requires an attached
// interpreter through `NativeStructureSource`.
unsafe impl Send for NativeSourceV2 {}

#[repr(C)]
pub(super) struct NativeCapsule {
    pub(super) api: NativeSourceV2,
    pub(super) structure: molframe::Structure,
    pub(super) topology: NativeTopology,
    pub(super) encoded_bcif: OnceLock<Result<Vec<u8>, ()>>,
    pub(super) queries: Mutex<crate::query_cache::QueryCache>,
    pub(super) pending: Mutex<crate::query_cache::PendingRows>,
}

unsafe impl Send for NativeCapsule {}
