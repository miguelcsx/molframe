//! DESRES Molecular Structure `SQLite` format.

mod error;
mod model;
mod reader;
mod schema;
mod writer;

pub use error::DmsError;
pub use model::{DmsBond, DmsCell, DmsFrame, DmsParticle, DmsSystem, DmsTopology, DmsVersion};
pub use reader::read_dms;
pub use writer::write_dms;

#[cfg(test)]
#[path = "dms_tests.rs"]
mod tests;
