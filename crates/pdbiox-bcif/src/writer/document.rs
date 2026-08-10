//! `BinaryCIF` document and structure writing.

use crate::container::{EncodedBlock, EncodedCategory, EncodedColumn, EncodedFile};
use crate::{EncodedData, encode_floats, encode_integers, encode_strings};
use num_traits::ToPrimitive;
use pdbiox_cif::{CifValue, Document};
use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::io::InputBuffer;
use pdbiox_core::structure::Structure;

/// Writes a shared CIF document as deterministic named-field `MessagePack`.
///
/// Each column tries the lossless candidates applicable to its value type and
/// keeps the smallest serialized representation. Equal-size ties are stable.
///
/// # Errors
///
/// Returns a registered codec diagnostic when a column is heterogeneous or
/// exceeds `BinaryCIF`'s integer representation.
pub fn write_document(document: &Document) -> Result<Vec<u8>, Diagnostic> {
    let mut data_blocks = Vec::with_capacity(document.len());
    for block in document.blocks() {
        let mut categories = Vec::with_capacity(block.len());
        for category in block.categories() {
            let mut columns = Vec::with_capacity(category.len());
            for item in category.items() {
                let Some(column) = category.column(item) else {
                    continue;
                };
                let (data, mask) = encode_column(column.iter())?;
                columns.push(EncodedColumn {
                    name: item.to_owned(),
                    data,
                    mask,
                });
            }
            categories.push(EncodedCategory {
                name: format!("_{}", category.name()),
                row_count: category.row_count(),
                columns,
            });
        }
        data_blocks.push(EncodedBlock {
            header: block.name().to_owned(),
            categories,
        });
    }
    let file = EncodedFile {
        version: "0.3.0".to_owned(),
        encoder: format!("pdbiox {}", env!("CARGO_PKG_VERSION")),
        data_blocks,
    };
    rmp_serde::to_vec_named(&file).map_err(|error| encode_error(&error))
}

/// Writes a structure by reusing the canonical mmCIF projection and parser.
///
/// This intentionally delegates structure semantics to `pdbiox-cif`; the
/// `BinaryCIF` crate owns only binary column encoding.
///
/// # Errors
///
/// Returns canonical projection, parsing, or `BinaryCIF` encoding diagnostics.
pub fn write_structure(structure: &Structure) -> Result<Vec<u8>, Vec<Diagnostic>> {
    write_structure_with_options(structure, &pdbiox_cif::CifWriteOptions::new())
}

/// Writes a structure with explicit canonical CIF identifier decisions.
///
/// # Errors
///
/// Returns canonical projection, parsing, or `BinaryCIF` encoding diagnostics.
pub fn write_structure_with_options(
    structure: &Structure,
    options: &pdbiox_cif::CifWriteOptions,
) -> Result<Vec<u8>, Vec<Diagnostic>> {
    let text = pdbiox_cif::write_canonical_with_options(structure, options).map_err(|error| {
        vec![
            Diagnostic::new(Code::E4105)
                .with_message("structure cannot be projected to canonical BinaryCIF")
                .with_context("reason", error.to_string()),
        ]
    })?;
    let input = InputBuffer::from_bytes(text.into_bytes());
    let (document, findings) = pdbiox_cif::parse(&input)?;
    if !findings.is_empty() {
        return Err(findings);
    }
    write_document(&document).map_err(|finding| vec![finding])
}

fn encode_column<'a>(
    values: impl Iterator<Item = &'a CifValue>,
) -> Result<(EncodedData, Option<EncodedData>), Diagnostic> {
    let values: Vec<&CifValue> = values.collect();
    let kind = column_kind(&values)?;
    let mut mask = Vec::with_capacity(values.len());
    let data = match kind {
        ColumnKind::Integer => {
            let mut output = Vec::with_capacity(values.len());
            for value in &values {
                output.push(integer_value(value, &mut mask));
            }
            encode_integers(&output)?
        }
        ColumnKind::Float => {
            let mut output = Vec::with_capacity(values.len());
            for value in &values {
                output.push(float_value(value, &mut mask)?);
            }
            encode_floats(&output)?
        }
        ColumnKind::String => {
            let mut output = Vec::with_capacity(values.len());
            for value in &values {
                output.push(string_value(value, &mut mask)?);
            }
            encode_strings(&output)?
        }
    };
    let has_mask = mask.iter().any(|value| *value != 0);
    let mask = if has_mask {
        Some(encode_integers(&mask)?)
    } else {
        None
    };
    Ok((data, mask))
}

#[derive(Clone, Copy)]
enum ColumnKind {
    Integer,
    Float,
    String,
}

fn column_kind(values: &[&CifValue]) -> Result<ColumnKind, Diagnostic> {
    let mut kind = None;
    for value in values {
        let candidate = match value {
            CifValue::Integer(_) => ColumnKind::Integer,
            CifValue::Float(_) => ColumnKind::Float,
            CifValue::Text(_) => ColumnKind::String,
            CifValue::Inapplicable | CifValue::Unknown => continue,
        };
        kind = match (kind, candidate) {
            (None, candidate) => Some(candidate),
            (Some(ColumnKind::Integer), ColumnKind::Integer) => Some(ColumnKind::Integer),
            (Some(ColumnKind::Integer), ColumnKind::Float)
            | (Some(ColumnKind::Float), ColumnKind::Integer | ColumnKind::Float) => {
                Some(ColumnKind::Float)
            }
            (Some(ColumnKind::String), ColumnKind::String) => Some(ColumnKind::String),
            _ => return Err(type_error()),
        };
    }
    Ok(match kind {
        Some(kind) => kind,
        None => ColumnKind::String,
    })
}

fn integer_value(value: &CifValue, mask: &mut Vec<i64>) -> i64 {
    match value {
        CifValue::Integer(value) => {
            mask.push(0);
            *value
        }
        CifValue::Inapplicable => {
            mask.push(1);
            0
        }
        _ => {
            mask.push(2);
            0
        }
    }
}

fn float_value(value: &CifValue, mask: &mut Vec<i64>) -> Result<f64, Diagnostic> {
    match value {
        CifValue::Integer(value) => {
            mask.push(0);
            value.to_f64().ok_or_else(type_error)
        }
        CifValue::Float(value) => {
            mask.push(0);
            Ok(*value)
        }
        CifValue::Inapplicable => {
            mask.push(1);
            Ok(0.0)
        }
        _ => {
            mask.push(2);
            Ok(0.0)
        }
    }
}

fn string_value(value: &CifValue, mask: &mut Vec<i64>) -> Result<String, Diagnostic> {
    match value {
        CifValue::Text(value) => {
            mask.push(0);
            Ok(value.to_string())
        }
        CifValue::Inapplicable => {
            mask.push(1);
            Ok(String::new())
        }
        CifValue::Unknown => {
            mask.push(2);
            Ok(String::new())
        }
        CifValue::Integer(_) | CifValue::Float(_) => Err(type_error()),
    }
}

fn type_error() -> Diagnostic {
    Diagnostic::new(Code::E1403).with_message("a BinaryCIF column mixes incompatible value types")
}

fn encode_error(error: &rmp_serde::encode::Error) -> Diagnostic {
    Diagnostic::new(Code::E1401).with_message(error.to_string())
}

#[cfg(test)]
#[path = "document_tests.rs"]
mod tests;
