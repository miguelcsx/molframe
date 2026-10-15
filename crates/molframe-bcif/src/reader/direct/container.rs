//! Borrowing container schema for the normal structure-read path.
//!
//! Large byte payloads point into the caller's input buffer. The public lazy
//! document remains fully owning, while a structure read avoids copying the
//! complete encoded file before column decoding begins.

use crate::codec::{
    Decoded, Encoding, Interner, check_declared_sizes, decode_borrowed, decode_shared, decoded_len,
    length_error,
};
use crate::container::{container_error, masked_value};
use molframe_cif::Category;
use molframe_cif::lexer::Quoting;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::io::Limits;
use molframe_core::span::ByteSpan;
use serde::Deserialize;

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
        for column in &category.columns {
            check_declared_sizes(&column.data.encoding, category.row_count)?;
            if let Some(mask) = &column.mask {
                check_declared_sizes(&mask.encoding, category.row_count)?;
            }
        }
    }
    Ok(file)
}

pub(super) fn decode_category(encoded: DirectCategory<'_>) -> Result<Category, Diagnostic> {
    let mut category = Category::new(encoded.name.trim_start_matches('_'), ByteSpan::default());
    let mut interner = Interner::default();
    for column in encoded.columns {
        let values = decode_shared(&column.data.encoding, column.data.data, &mut interner)?;
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

#[cfg(test)]
#[path = "container_tests.rs"]
mod tests;
