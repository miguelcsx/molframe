//! `BinaryCIF` codec types and decoding.

mod bounds;
mod decode;
mod numeric_into;
mod string_column;

pub(crate) use bounds::check_declared_sizes;
pub(crate) use decode::decode_borrowed;
pub use decode::*;
pub(crate) use numeric_into::decode_f32_borrowed_into;
pub use numeric_into::decode_f32_into;
pub use string_column::{DecodedStringColumn, DecodedStringIter};
