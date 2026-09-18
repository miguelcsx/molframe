//! Shared deterministic encoded-size selection.

use crate::codec::{EncodedData, Encoding};

pub(super) fn encoded_size(candidate: &EncodedData) -> usize {
    match crate::messagepack::compact_size(candidate) {
        Ok(bytes) => bytes,
        Err(_) => usize::MAX,
    }
}

pub(super) fn encoded_size_lower_bound(encoding: Vec<Encoding>, data_len: usize) -> Option<usize> {
    let empty = EncodedData {
        encoding,
        data: Vec::new(),
    };
    let empty_size = crate::messagepack::compact_size(&empty).ok()?;
    empty_size
        .checked_sub(byte_string_size(0)?)?
        .checked_add(byte_string_size(data_len)?)
}

pub(super) struct EncodingChoice {
    value: EncodedData,
    size: usize,
}

impl EncodingChoice {
    pub(super) fn new(value: EncodedData) -> Self {
        let size = encoded_size(&value);
        Self { value, size }
    }

    pub(super) const fn size(&self) -> usize {
        self.size
    }

    pub(super) fn consider(&mut self, candidate: EncodedData) {
        let size = encoded_size(&candidate);
        if size < self.size {
            self.value = candidate;
            self.size = size;
        }
    }

    pub(super) fn finish(self) -> EncodedData {
        self.value
    }
}

fn byte_string_size(length: usize) -> Option<usize> {
    let header = if u8::try_from(length).is_ok() {
        2
    } else if u16::try_from(length).is_ok() {
        3
    } else if length <= usize::try_from(u32::MAX).ok()? {
        5
    } else {
        return None;
    };
    length.checked_add(header)
}

#[cfg(test)]
#[path = "strategy_tests.rs"]
mod tests;
