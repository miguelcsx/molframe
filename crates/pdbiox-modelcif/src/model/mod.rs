//! Typed `ModelCIF` data model and structure extension.

mod category;
mod confidence;
mod error;
mod metadata;
mod packed;
mod schema;
mod view;

pub use category::*;
pub use confidence::*;
pub use error::*;
pub use metadata::*;
pub(crate) use packed::{PackedIndices, PackedIntegers};
pub(crate) use schema::identifier_item;
pub use view::*;
