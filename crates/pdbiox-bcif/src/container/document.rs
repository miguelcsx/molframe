//! `MessagePack` container and lazy category decoding.

use crate::codec::{Decoded, EncodedData, decode};
use pdbiox_cif::lexer::Quoting;
use pdbiox_cif::{Category, CifValue, DataBlock, Document};
use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::io::Limits;
use pdbiox_core::span::ByteSpan;
use serde::{Deserialize, Serialize};
use std::sync::Arc;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub(crate) struct EncodedFile {
    pub(crate) version: String,
    pub(crate) encoder: String,
    #[serde(rename = "dataBlocks")]
    pub(crate) data_blocks: Vec<EncodedBlock>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub(crate) struct EncodedBlock {
    pub(crate) header: String,
    pub(crate) categories: Vec<EncodedCategory>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub(crate) struct EncodedCategory {
    pub(crate) name: String,
    #[serde(rename = "rowCount")]
    pub(crate) row_count: usize,
    pub(crate) columns: Vec<EncodedColumn>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub(crate) struct EncodedColumn {
    pub(crate) name: String,
    pub(crate) data: EncodedData,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) mask: Option<EncodedData>,
}

/// A parsed `BinaryCIF` container whose columns remain encoded until requested.
#[derive(Clone, Debug)]
pub struct BinaryDocument {
    pub(crate) file: EncodedFile,
    limits: Limits,
}

impl BinaryDocument {
    /// Parses only the `MessagePack` container and encoding metadata.
    ///
    /// # Errors
    ///
    /// Returns a registered encoding or length diagnostic for malformed input.
    pub fn parse(bytes: &[u8], limits: Limits) -> Result<Self, Diagnostic> {
        let byte_count =
            u64::try_from(bytes.len()).map_err(|_| Limits::exceeded("input bytes", bytes.len()))?;
        if byte_count > limits.decompressed_bytes {
            return Err(Limits::exceeded("input bytes", bytes.len()));
        }
        let file: EncodedFile =
            rmp_serde::from_slice(bytes).map_err(|error| container_error(&error))?;
        for category in file.data_blocks.iter().flat_map(|block| &block.categories) {
            let row_count = u64::try_from(category.row_count)
                .map_err(|_| Limits::exceeded("rows per category", category.row_count))?;
            if row_count > limits.rows_per_category {
                return Err(Limits::exceeded("rows per category", category.row_count));
            }
        }
        Ok(Self { file, limits })
    }

    /// Format version declared by the container.
    #[must_use]
    pub fn version(&self) -> &str {
        &self.file.version
    }

    /// Encoder identity declared by the container.
    #[must_use]
    pub fn encoder(&self) -> &str {
        &self.file.encoder
    }

    /// Number of data blocks without decoding any column.
    #[must_use]
    pub fn block_count(&self) -> usize {
        self.file.data_blocks.len()
    }

    /// Decodes one named category while leaving every other category encoded.
    ///
    /// # Errors
    ///
    /// Returns codec diagnostics and declared-length mismatches.
    pub fn category(&self, block: usize, name: &str) -> Result<Option<Category>, Diagnostic> {
        let Some(block) = self.file.data_blocks.get(block) else {
            return Ok(None);
        };
        let name = name.trim_start_matches('_');
        let Some(category) = block
            .categories
            .iter()
            .find(|category| category.name.trim_start_matches('_') == name)
        else {
            return Ok(None);
        };
        decode_category(category).map(Some)
    }

    /// Decodes every category into the shared lossless document model.
    ///
    /// # Errors
    ///
    /// Returns the first malformed column in stable block/category order.
    pub fn to_document(&self) -> Result<Document, Diagnostic> {
        let mut document = Document::new();
        for encoded_block in &self.file.data_blocks {
            let mut block = DataBlock::new(encoded_block.header.clone());
            for encoded_category in &encoded_block.categories {
                let decoded = decode_category(encoded_category)?;
                let name = decoded.name().to_owned();
                *block.category_mut(&name, ByteSpan::default()) = decoded;
            }
            document.push(block);
        }
        Ok(document)
    }

    /// Limits retained for subsequent lazy decoding.
    #[must_use]
    pub const fn limits(&self) -> Limits {
        self.limits
    }
}

fn decode_category(encoded: &EncodedCategory) -> Result<Category, Diagnostic> {
    let mut category = Category::new(encoded.name.trim_start_matches('_'), ByteSpan::default());
    for column in &encoded.columns {
        let values = decode(&column.data)?;
        if decoded_len(&values) != encoded.row_count {
            return Err(length_error(decoded_len(&values), encoded.row_count));
        }
        let mask = match &column.mask {
            Some(mask) => {
                let Decoded::Integers(mask) = decode(mask)? else {
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
        if let Some(mask) = mask {
            for row in 0..encoded.row_count {
                let value = masked_value(&values, Some(&mask), row)?;
                target.push(value, Quoting::Bare);
            }
        } else {
            match values {
                Decoded::Integers(values) => {
                    for value in values {
                        target.push(CifValue::Integer(value), Quoting::Bare);
                    }
                }
                Decoded::Floats(values) => {
                    for value in values {
                        target.push(CifValue::Float(value), Quoting::Bare);
                    }
                }
                Decoded::Strings(values) => {
                    for row in 0..values.len() {
                        let Some(value) = values.get_shared(row) else {
                            return Err(length_error(row, values.len()));
                        };
                        target.push(CifValue::Text(Arc::clone(value)), Quoting::Bare);
                    }
                }
            }
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
