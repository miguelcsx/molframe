//! Progressive multiple-sequence alignment with optional iterative refinement.

mod algorithm;
mod guide;
mod profile;
mod types;

pub use algorithm::{msa, progressive_msa};
pub use types::{MsaError, MsaOptions};
