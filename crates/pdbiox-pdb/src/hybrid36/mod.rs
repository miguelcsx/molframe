//! Hybrid-36 fixed-width integer coding.

mod codec;

mod compatibility {
    pub use super::codec::needs_encoding as needs_hybrid36;
}

pub use codec::*;
pub use compatibility::*;
