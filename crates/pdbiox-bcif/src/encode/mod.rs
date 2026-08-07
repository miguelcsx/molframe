//! Deterministic `BinaryCIF` encoders and per-column optimisation.

mod float;
mod integer;
mod string;

pub use float::{encode_floats, encode_interval};
pub use integer::encode_integers;
pub use string::encode_strings;

use crate::codec::EncodedData;

fn encoded_size(candidate: &EncodedData) -> usize {
    match rmp_serde::to_vec(candidate) {
        Ok(bytes) => bytes.len(),
        Err(_) => usize::MAX,
    }
}

fn keep_smaller(current: &mut EncodedData, candidate: EncodedData) {
    if encoded_size(&candidate) < encoded_size(current) {
        *current = candidate;
    }
}
