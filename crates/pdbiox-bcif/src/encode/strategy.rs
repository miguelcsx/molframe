//! Shared deterministic encoded-size selection.

use crate::codec::EncodedData;

fn encoded_size(candidate: &EncodedData) -> usize {
    match rmp_serde::to_vec(candidate) {
        Ok(bytes) => bytes.len(),
        Err(_) => usize::MAX,
    }
}

pub(super) fn keep_smaller(current: &mut EncodedData, candidate: EncodedData) {
    if encoded_size(&candidate) < encoded_size(current) {
        *current = candidate;
    }
}
