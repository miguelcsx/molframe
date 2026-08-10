//! Deterministic `BinaryCIF` encoders and per-column optimisation.

mod float;
mod integer;
mod strategy;
mod string;

pub use float::{encode_floats, encode_interval};
pub use integer::encode_integers;
pub use string::encode_strings;
