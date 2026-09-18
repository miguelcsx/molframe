//! Complete formatted AMBER restart records and writing.

mod model;
mod write;

pub use model::AmberRestart;
pub use write::write_amber_restart;
