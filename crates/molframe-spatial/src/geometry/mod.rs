//! Periodic geometry and numeric boundaries.

pub(crate) mod numeric;
pub(crate) mod periodic;

pub use periodic::{PeriodicBox, PeriodicImage};
