//! `BinaryCIF` codec types and decoding.

mod bounds;
mod decode;
mod interner;
mod numeric_into;
mod string_column;
mod strings;

pub(crate) use bounds::check_declared_sizes;
pub use decode::{DataType, Decoded, EncodedData, Encoding, decode};
pub(crate) use decode::{decode_borrowed, decode_shared, decoded_len, length_error};
pub(crate) use interner::Interner;
pub(crate) use numeric_into::decode_f32_borrowed_into;
pub use numeric_into::decode_f32_into;
pub use string_column::{DecodedStringColumn, DecodedStringIter};
