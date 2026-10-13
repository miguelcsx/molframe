//! Dense confidence storage: budgeted PAE matrices and their precision layout.

// The storage codec converts between widths by design; every cast is bounded
// by the layout chosen from the observed values and covered by round-trip tests.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::float_cmp
)]

mod error;
mod half;
mod json;
mod layout;
mod pae;
mod precision;

pub use error::ConfidenceError;
pub use json::read_pae_json;
pub use layout::{PaeKind, PaeLayout, PaePolicy, PaeStorage};
pub use pae::PaeMatrix;
pub use precision::PrecisionProbe;
