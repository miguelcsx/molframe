//! Hybrid-36 fixed-width integer coding.

mod codec;

pub use codec::needs_encoding as needs_hybrid36;
pub use codec::{decode, encode, needs_encoding};
