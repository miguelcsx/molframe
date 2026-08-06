//! Column storage: encodings, bit-level primitives, and why a value is absent.

mod bits;
mod encoded;
mod validity;

pub use bits::{BitVec, bit_width, pack, unpack_one};
pub use encoded::{ColumnIter, ColumnValue, EncodedColumn};
pub use validity::{Presence, ValidityMask};
