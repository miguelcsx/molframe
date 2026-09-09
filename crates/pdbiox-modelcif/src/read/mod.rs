//! Direct allocation-bounded `ModelCIF` projection from CIF parser events.

mod batch;
mod category;
mod column;
mod packed;
mod projection;
mod reader;

pub use batch::ModelCifBatchSource;
pub use projection::ModelCifProjection;
pub use reader::{read_compact, read_compact_with_options};

#[cfg(test)]
#[path = "read_tests.rs"]
mod tests;
