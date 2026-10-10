//! Periodic geometry and numeric boundaries.

pub(crate) mod numeric;
pub(crate) mod periodic;
mod reduction;

pub use periodic::{PeriodicBox, PeriodicImage};
