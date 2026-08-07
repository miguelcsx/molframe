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
    let mut dictionary = IndexMap::<&str, i64>::new();
    let mut indices = Vec::with_capacity(values.len());
    for value in values {
        let index = if let Some(index) = dictionary.get(value.as_str()) {
            *index
        } else {
            let Ok(index) = i64::try_from(dictionary.len()) else {
                return Err(length_error());
            };
            dictionary.insert(value, index);
            index
        };
        indices.push(index);
    }

    let mut string_data = String::new();
    let mut offsets = Vec::with_capacity(dictionary.len().saturating_add(1));
    offsets.push(0);
    for value in dictionary.keys() {
        string_data.push_str(value);
        let Ok(offset) = i64::try_from(string_data.len()) else {
            return Err(length_error());
        };
        offsets.push(offset);
    }

    let indices = encode_integers(&indices)?;
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
