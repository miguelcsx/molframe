//! Direct `BinaryCIF` structure reading without a coordinate DOM.
//!
//! Container metadata and the columns consumed by structure lowering are
//! decoded once. Coordinate rows then borrow dictionary strings and numeric
//! arrays while feeding the shared CIF lowerer synchronously.

mod column;
mod container;
mod ensemble;
mod external;
mod projection;
mod reader;
mod values;

pub(super) use reader::{read, read_with_metadata, read_with_projection};
