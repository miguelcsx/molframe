//! `BinaryCIF` document and structure reading.

mod batch;
mod direct;
mod index;
mod structure;

pub use batch::BcifBatchSource;
pub use structure::*;
