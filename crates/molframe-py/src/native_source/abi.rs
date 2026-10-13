//! The C layout shared by a structure's producer and its consumers.

use std::ffi::{CStr, c_char};
use std::sync::{Mutex, OnceLock};

pub(super) const CAPSULE_NAME: &CStr = c"molframe.StructureSource.v2";
/// The version this producer writes: the base layout plus the failure channel.
pub(super) const ABI_VERSION: u32 = 3;
/// The oldest layout a consumer still reads; it has no failure channel.
pub(super) const BASE_ABI_VERSION: u32 = 2;

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
    /// Compact `MolFrame` bond-order code.
    ///
    /// The byte occupied the first slot of the old reserved tail, so the
    /// capsule layout remains compatible with consumers that only read
    /// endpoints and aromaticity.
    pub order: u8,
    /// One only when the source explicitly records metal coordination.
    /// This occupies a formerly reserved byte without changing the C layout.
    pub metal_coordination: u8,
    pub(super) reserved: [u8; 1],
}

impl NativeBond {
    /// Unknown or unavailable source order.
    pub const ORDER_UNKNOWN: u8 = 0;
    /// Single source order.
    pub const ORDER_SINGLE: u8 = 1;
    /// Double source order.
    pub const ORDER_DOUBLE: u8 = 2;
    /// Triple source order.
    pub const ORDER_TRIPLE: u8 = 3;
    /// Quadruple source order.
    pub const ORDER_QUADRUPLE: u8 = 4;
    /// Aromatic source order.
    pub const ORDER_AROMATIC: u8 = 5;
    /// Polymeric source order.
    pub const ORDER_POLYMERIC: u8 = 6;
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
/// Reports why the last BCIF encoding failed: `(source, code, message, capacity)`.
///
/// Writes the diagnostic number (zero when the failure names none) through
/// `code` when it is non-null, and returns the UTF-8 message length. A null
/// `message` or zero `capacity` asks only for the length. Returns a negative
/// value when no failure is recorded or the buffer is too small.
pub(super) type FailureFn =
    unsafe extern "C" fn(*const NativeSourceV2, *mut u32, *mut u8, usize) -> i64;

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

/// The version-3 layout: the base record followed by the failure channel.
///
/// Appending leaves every earlier field where a version-2 consumer expects
/// it; a consumer reads `describe_failure` only after the version check.
#[repr(C)]
pub(super) struct NativeSourceV3 {
    pub(super) base: NativeSourceV2,
    pub(super) describe_failure: FailureFn,
}

/// Why an encoding failed, kept beside the cached result.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct NativeFailure {
    /// The `MolFrame` diagnostic number, or zero when none was named.
    pub(super) code: u32,
    /// The diagnostic's own message.
    pub(super) message: String,
}

// The raw pointers refer only to immutable buffers retained by the same
// capsule. Python owns capsule destruction, and access requires an attached
// interpreter through `NativeStructureSource`.
unsafe impl Send for NativeSourceV2 {}

#[repr(C)]
pub(super) struct NativeCapsule {
    pub(super) api: NativeSourceV3,
    pub(super) structure: molframe::Structure,
    pub(super) topology: NativeTopology,
    pub(super) encoded_bcif: OnceLock<Result<Vec<u8>, NativeFailure>>,
    pub(super) queries: Mutex<crate::query_cache::QueryCache>,
    pub(super) pending: Mutex<crate::query_cache::PendingRows>,
}

unsafe impl Send for NativeCapsule {}
