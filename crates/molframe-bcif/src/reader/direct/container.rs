//! Borrowing container schema for the normal structure-read path.
//!
//! Large byte payloads point into the caller's input buffer. The public lazy
//! document remains fully owning, while a structure read avoids copying the
//! complete encoded file before column decoding begins.

use crate::codec::{Decoded, Encoding, decode_borrowed};
use molframe_cif::lexer::Quoting;
use molframe_cif::{Category, CifValue};
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::io::Limits;
use molframe_core::span::ByteSpan;
use serde::Deserialize;
use std::sync::Arc;

#[derive(Deserialize)]
pub(super) struct DirectFile<'a> {
    #[serde(rename = "dataBlocks", borrow)]
    pub(super) data_blocks: Vec<DirectBlock<'a>>,
}

#[derive(Deserialize)]
pub(super) struct DirectBlock<'a> {
    pub(super) header: String,
    #[serde(borrow)]
    pub(super) categories: Vec<DirectCategory<'a>>,
}

#[derive(Deserialize)]
pub(super) struct DirectCategory<'a> {
    pub(super) name: String,
    #[serde(rename = "rowCount")]
    pub(super) row_count: usize,
    #[serde(borrow)]
    pub(super) columns: Vec<DirectColumn<'a>>,
}

#[derive(Deserialize)]
pub(super) struct DirectColumn<'a> {
    pub(super) name: String,
    #[serde(borrow)]
    pub(super) data: DirectData<'a>,
    #[serde(default, borrow)]
    pub(super) mask: Option<DirectData<'a>>,
}

#[derive(Deserialize)]
pub(super) struct DirectData<'a> {
    pub(super) encoding: Vec<Encoding>,
    #[serde(borrow, with = "serde_bytes")]
    pub(super) data: &'a [u8],
}

pub(super) fn parse(bytes: &[u8], limits: Limits) -> Result<DirectFile<'_>, Diagnostic> {
    let byte_count =
        u64::try_from(bytes.len()).map_err(|_| Limits::exceeded("input bytes", bytes.len()))?;
    if byte_count > limits.decompressed_bytes {
        return Err(Limits::exceeded("input bytes", bytes.len()));
    }
    let file: DirectFile<'_> =
        rmp_serde::from_slice(bytes).map_err(|error| container_error(&error))?;
    for category in file.data_blocks.iter().flat_map(|block| &block.categories) {
        let row_count = u64::try_from(category.row_count)
            .map_err(|_| Limits::exceeded("rows per category", category.row_count))?;
        if row_count > limits.rows_per_category {
            return Err(Limits::exceeded("rows per category", category.row_count));
        }
    }
    Ok(file)
}

pub(super) fn decode_category(encoded: DirectCategory<'_>) -> Result<Category, Diagnostic> {
    let mut category = Category::new(encoded.name.trim_start_matches('_'), ByteSpan::default());
    for column in encoded.columns {
        let values = decode_borrowed(&column.data.encoding, column.data.data)?;
        if decoded_len(&values) != encoded.row_count {
            return Err(length_error(decoded_len(&values), encoded.row_count));
        }
        let mask = match column.mask {
            Some(mask) => {
                let Decoded::Integers(mask) = decode_borrowed(&mask.encoding, mask.data)? else {
                    return Err(Diagnostic::new(Code::E1403)
                        .with_message("a BinaryCIF mask is not an integer array"));
                };
                if mask.len() != encoded.row_count {
                    return Err(length_error(mask.len(), encoded.row_count));
                }
                Some(mask)
            }
            None => None,
        };
        let target = category.column_mut(&column.name);
        target.reserve(encoded.row_count);
        for row in 0..encoded.row_count {
            target.push(masked_value(&values, mask.as_deref(), row)?, Quoting::Bare);
        }
    }
    Ok(category)
}

fn masked_value(
    values: &Decoded,
    mask: Option<&[i64]>,
    row: usize,
) -> Result<CifValue, Diagnostic> {
    match mask.and_then(|values| values.get(row)).copied() {
        Some(1) => return Ok(CifValue::Inapplicable),
        Some(2) => return Ok(CifValue::Unknown),
        Some(0) | None => {}
        Some(value) => {
            return Err(Diagnostic::new(Code::E1401).with_context("mask value", value.to_string()));
        }
    }
    match values {
        Decoded::Integers(values) => Ok(CifValue::Integer(
            *values
                .get(row)
                .ok_or_else(|| length_error(row, values.len()))?,
        )),
        Decoded::Floats(values) => Ok(CifValue::Float(
            *values
                .get(row)
                .ok_or_else(|| length_error(row, values.len()))?,
        )),
        Decoded::Strings(values) => Ok(CifValue::Text(Arc::clone(
            values
                .get_shared(row)
                .ok_or_else(|| length_error(row, values.len()))?,
        ))),
    }
}

fn decoded_len(values: &Decoded) -> usize {
    match values {
        Decoded::Integers(values) => values.len(),
        Decoded::Floats(values) => values.len(),
        Decoded::Strings(values) => values.len(),
    }
}

fn container_error(error: &rmp_serde::decode::Error) -> Diagnostic {
    let message = error.to_string();
    let code = if message.contains("unknown variant") {
        Code::E1402
    } else {
        Code::E1401
    };
    Diagnostic::new(code).with_message(message)
}

fn length_error(actual: usize, declared: usize) -> Diagnostic {
    Diagnostic::new(Code::E1401)
        .with_context("actual", actual.to_string())
        .with_context("declared", declared.to_string())
}

#[cfg(test)]
#[path = "container_tests.rs"]
mod tests;
