//! Stable dictionary encoding for text columns.

use super::integer::encode_integers;
use crate::codec::{EncodedData, Encoding};
use indexmap::IndexMap;
use pdbiox_core::diagnostic::{Code, Diagnostic};

/// Encodes strings using first-seen dictionary order and encoded offsets.
///
/// # Errors
///
/// Returns a length diagnostic when dictionary offsets exceed `BinaryCIF`'s
/// integer representation.
pub fn encode_strings(values: &[String]) -> Result<EncodedData, Diagnostic> {
    let mut encoder = StringColumnEncoder::with_capacity(values.len());
    for value in values {
        encoder.push(value);
    }
    encoder.finish()
}

pub(crate) struct StringColumnEncoder {
    dictionary: IndexMap<Box<str>, i64>,
    indices: Vec<i64>,
    invalid: bool,
}

impl StringColumnEncoder {
    pub(crate) fn with_capacity(rows: usize) -> Self {
        Self {
            dictionary: IndexMap::new(),
            indices: Vec::with_capacity(rows),
            invalid: false,
        }
    }

    pub(crate) fn push(&mut self, value: &str) {
        let index = if let Some(index) = self.dictionary.get(value) {
            *index
        } else {
            let Ok(index) = i64::try_from(self.dictionary.len()) else {
                self.invalid = true;
                self.indices.push(0);
                return;
            };
            self.dictionary.insert(value.into(), index);
            index
        };
        self.indices.push(index);
    }

    pub(crate) fn push_ascii_uppercase(&mut self, value: &str) {
        let bytes = value.as_bytes();
        let mut uppercase = [0_u8; 3];
        if bytes.len() > uppercase.len() || bytes.iter().any(|byte| !byte.is_ascii()) {
            self.invalid = true;
            self.indices.push(0);
            return;
        }
        for (target, source) in uppercase.iter_mut().zip(bytes) {
            *target = source.to_ascii_uppercase();
        }
        let Ok(text) = std::str::from_utf8(&uppercase[..bytes.len()]) else {
            self.invalid = true;
            self.indices.push(0);
            return;
        };
        self.push(text);
    }

    pub(crate) fn finish(self) -> Result<EncodedData, Diagnostic> {
        if self.invalid {
            return Err(length_error());
        }
        encode_dictionary(&self.dictionary, &self.indices)
    }
}

fn encode_dictionary(
    dictionary: &IndexMap<Box<str>, i64>,
    indices: &[i64],
) -> Result<EncodedData, Diagnostic> {
    let string_capacity = dictionary.keys().try_fold(0usize, |total, value| {
        total.checked_add(value.len()).ok_or_else(length_error)
    })?;
    let mut string_data = String::with_capacity(string_capacity);
    let Some(offset_capacity) = dictionary.len().checked_add(1) else {
        return Err(length_error());
    };
    let mut offsets = Vec::with_capacity(offset_capacity);
    offsets.push(0);
    for value in dictionary.keys() {
        string_data.push_str(value.as_ref());
        let Ok(offset) = i64::try_from(string_data.len()) else {
            return Err(length_error());
        };
        offsets.push(offset);
    }

    let indices = encode_integers(indices)?;
    let offsets = encode_integers(&offsets)?;
    Ok(EncodedData {
        encoding: vec![Encoding::StringArray {
            data_encoding: indices.encoding,
            string_data,
            offset_encoding: offsets.encoding,
            offsets: offsets.data,
        }],
        data: indices.data,
    })
}

fn length_error() -> Diagnostic {
    Diagnostic::new(Code::E1401).with_message("string dictionary is too large")
}

#[cfg(test)]
#[path = "string_tests.rs"]
mod tests;
