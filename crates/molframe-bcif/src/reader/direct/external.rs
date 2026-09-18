//! Column-at-a-time delivery to an external compact projection.

use super::container::{DirectCategory, DirectData};
use crate::codec::{Decoded, decode_borrowed};
use molframe_cif::{Category, CifEventSink, CifScalar, CifValue};
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::span::ByteSpan;

pub(super) fn feed_encoded<S>(encoded: DirectCategory<'_>, sink: &mut S) -> Result<(), Diagnostic>
where
    S: CifEventSink,
{
    let DirectCategory {
        name,
        row_count,
        columns,
    } = encoded;
    let category = name.trim_start_matches('_');
    let tags = columns
        .iter()
        .map(|column| (category.into(), column.name.as_str().into()))
        .collect::<Vec<_>>();
    sink.begin_loop(&tags);
    for column in columns {
        if !sink.accepts_category(category) {
            break;
        }
        let values = decode_borrowed(&column.data.encoding, column.data.data)?;
        if decoded_len(&values) != row_count {
            return Err(length_error(decoded_len(&values), row_count));
        }
        let mask = decode_mask(column.mask, row_count)?;
        for row in 0..row_count {
            sink.value(
                category,
                &column.name,
                scalar(&values, mask.as_deref(), row)?,
                ByteSpan::default(),
            );
        }
    }
    for _ in 0..row_count {
        sink.end_row();
    }
    sink.end_loop();
    Ok(())
}

pub(super) fn feed_decoded<S>(category: &Category, sink: &mut S)
where
    S: CifEventSink,
{
    let tags = category
        .items()
        .map(|item| (category.name().into(), item.into()))
        .collect::<Vec<_>>();
    sink.begin_loop(&tags);
    for item in category.items() {
        for row in 0..category.row_count() {
            let value = match category.value(item, row) {
                Some(value) => scalar_from_value(value),
                None => CifScalar::Unknown,
            };
            sink.value(category.name(), item, value, ByteSpan::default());
        }
    }
    for _ in 0..category.row_count() {
        sink.end_row();
    }
    sink.end_loop();
}

fn decode_mask(
    encoded: Option<DirectData<'_>>,
    row_count: usize,
) -> Result<Option<Vec<i64>>, Diagnostic> {
    let Some(encoded) = encoded else {
        return Ok(None);
    };
    let Decoded::Integers(mask) = decode_borrowed(&encoded.encoding, encoded.data)? else {
        return Err(
            Diagnostic::new(Code::E1403).with_message("a BinaryCIF mask is not an integer array")
        );
    };
    if mask.len() != row_count {
        return Err(length_error(mask.len(), row_count));
    }
    Ok(Some(mask))
}

fn scalar<'a>(
    values: &'a Decoded,
    mask: Option<&[i64]>,
    row: usize,
) -> Result<CifScalar<'a>, Diagnostic> {
    match mask.and_then(|values| values.get(row)).copied() {
        Some(1) => return Ok(CifScalar::Inapplicable),
        Some(2) => return Ok(CifScalar::Unknown),
        Some(0) | None => {}
        Some(value) => {
            return Err(Diagnostic::new(Code::E1401).with_context("mask value", value.to_string()));
        }
    }
    match values {
        Decoded::Integers(values) => values
            .get(row)
            .copied()
            .map(CifScalar::Integer)
            .ok_or_else(|| length_error(row, values.len())),
        Decoded::Floats(values) => values
            .get(row)
            .copied()
            .map(CifScalar::Float)
            .ok_or_else(|| length_error(row, values.len())),
        Decoded::Strings(values) => values
            .get(row)
            .map(CifScalar::Text)
            .ok_or_else(|| length_error(row, values.len())),
    }
}

fn scalar_from_value(value: &CifValue) -> CifScalar<'_> {
    match value {
        CifValue::Inapplicable => CifScalar::Inapplicable,
        CifValue::Unknown => CifScalar::Unknown,
        CifValue::Text(value) => CifScalar::Text(value),
        CifValue::Integer(value) => CifScalar::Integer(*value),
        CifValue::Float(value) => CifScalar::Float(*value),
    }
}

fn decoded_len(values: &Decoded) -> usize {
    match values {
        Decoded::Integers(values) => values.len(),
        Decoded::Floats(values) => values.len(),
        Decoded::Strings(values) => values.len(),
    }
}

fn length_error(actual: usize, declared: usize) -> Diagnostic {
    Diagnostic::new(Code::E1401)
        .with_context("actual", actual.to_string())
        .with_context("declared", declared.to_string())
}
