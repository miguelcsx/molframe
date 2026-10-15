//! Trajectory Next Generation container support.

mod error;
mod model;
mod reader;
mod writer;

pub use error::TngError;
use model::length_to_angstrom;
pub use model::{TngCompression, TngTrajectory, TngWriteOptions};
pub use reader::parse_tng;
pub use writer::write_tng;

#[cfg(test)]
#[path = "tng_tests.rs"]
mod tests;
