//! Versioned native structure sharing between independently loaded extensions.
//!
//! The capsule contains only C-layout records, primitive slices and one C-ABI
//! query callback. Rust-owned values never cross the extension boundary. The
//! consumer retains the capsule, which in turn retains the immutable structure
//! snapshot and every exported buffer.
//!
//! `abi` holds the C layout both sides agree on, `consumer` is the safe view a
//! second extension imports, and `producer` builds the capsule and answers its
//! callbacks.

mod abi;
mod consumer;
mod producer;

#[cfg(test)]
#[path = "../native_source_tests.rs"]
mod tests;

pub use abi::{NativeAtom, NativeBond, NativeTopology};
pub use consumer::NativeStructureSource;
pub(crate) use producer::capsule;
