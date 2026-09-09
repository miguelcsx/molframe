//! Bounded `BinaryCIF` container index over stable source offsets.

mod cursor;
mod marker;

use cursor::{Cursor, error};
use marker::Value;
use pdbiox_core::{Diagnostic, ExecutionContext, MemoryReservation, SourceBytes};
use std::ops::Range;

const INDEX_BASE_BYTES: usize = 64 * 1024;
const INDEX_ENTRY_BYTES: usize = 512;

#[derive(Clone, Debug)]
pub(crate) struct IndexedData {
    pub(crate) encoding: Range<u64>,
    pub(crate) payload: Range<u64>,
}

#[derive(Clone, Debug)]
pub(crate) struct IndexedColumn {
    pub(crate) name: String,
    pub(crate) data: IndexedData,
    pub(crate) mask: Option<IndexedData>,
}

#[derive(Clone, Debug)]
pub(crate) struct IndexedCategory {
    pub(crate) name: String,
    pub(crate) rows: u64,
    pub(crate) columns: Vec<IndexedColumn>,
}

#[derive(Debug)]
pub(crate) struct BcifOffsetIndex {
    pub(crate) categories: Vec<IndexedCategory>,
    _reservation: MemoryReservation,
}

impl BcifOffsetIndex {
    pub(crate) fn build<S: SourceBytes>(
        source: &mut S,
        window_bytes: usize,
        context: &ExecutionContext,
    ) -> Result<Self, Diagnostic> {
        let mut reservation = context
            .try_reserve(INDEX_BASE_BYTES)
            .map_err(index_memory)?;
        let mut cursor = Cursor::new(source, window_bytes, context)?;
        let root = expect_map(marker::value(&mut cursor)?)?;
        let mut categories = Vec::new();
        for _ in 0..root {
            let key = read_key(&mut cursor)?;
            let value = marker::value(&mut cursor)?;
            if key == "dataBlocks" {
                parse_blocks(&mut cursor, value, &mut categories, &mut reservation)?;
            } else {
                marker::skip(&mut cursor, value)?;
            }
        }
        Ok(Self {
            categories,
            _reservation: reservation,
        })
    }
}

fn parse_blocks<S: SourceBytes>(
    cursor: &mut Cursor<'_, S>,
    value: Value,
    categories: &mut Vec<IndexedCategory>,
    reservation: &mut MemoryReservation,
) -> Result<(), Diagnostic> {
    let blocks = expect_array(value)?;
    for _ in 0..blocks {
        let fields = expect_map(marker::value(cursor)?)?;
        for _ in 0..fields {
            let key = read_key(cursor)?;
            let value = marker::value(cursor)?;
            if key == "categories" {
                parse_categories(cursor, value, categories, reservation)?;
            } else {
                marker::skip(cursor, value)?;
            }
        }
    }
    Ok(())
}

fn parse_categories<S: SourceBytes>(
    cursor: &mut Cursor<'_, S>,
    value: Value,
    categories: &mut Vec<IndexedCategory>,
    reservation: &mut MemoryReservation,
) -> Result<(), Diagnostic> {
    let count = expect_array(value)?;
    for _ in 0..count {
        reservation
            .try_grow(INDEX_ENTRY_BYTES)
            .map_err(index_memory)?;
        let fields = expect_map(marker::value(cursor)?)?;
        let mut category = IndexedCategory {
            name: String::new(),
            rows: 0,
            columns: Vec::new(),
        };
        for _ in 0..fields {
            let key = read_key(cursor)?;
            let value = marker::value(cursor)?;
            match key.as_str() {
                "name" => category.name = read_text(cursor, value)?,
                "rowCount" => category.rows = expect_unsigned(value)?,
                "columns" => {
                    parse_columns(cursor, value, &mut category.columns, reservation)?;
                }
                _ => marker::skip(cursor, value)?,
            }
        }
        categories.push(category);
    }
    Ok(())
}

fn parse_columns<S: SourceBytes>(
    cursor: &mut Cursor<'_, S>,
    value: Value,
    columns: &mut Vec<IndexedColumn>,
    reservation: &mut MemoryReservation,
) -> Result<(), Diagnostic> {
    let count = expect_array(value)?;
    for _ in 0..count {
        reservation
            .try_grow(INDEX_ENTRY_BYTES)
            .map_err(index_memory)?;
        let fields = expect_map(marker::value(cursor)?)?;
        let mut name = String::new();
        let mut data = None;
        let mut mask = None;
        for _ in 0..fields {
            let key = read_key(cursor)?;
            let value = marker::value(cursor)?;
            match key.as_str() {
                "name" => name = read_text(cursor, value)?,
                "data" => data = Some(parse_data(cursor, value)?),
                "mask" if matches!(value, Value::Nil) => {}
                "mask" => mask = Some(parse_data(cursor, value)?),
                _ => marker::skip(cursor, value)?,
            }
        }
        let Some(data) = data else {
            return Err(error("BinaryCIF column has no data object"));
        };
        columns.push(IndexedColumn { name, data, mask });
    }
    Ok(())
}

fn parse_data<S: SourceBytes>(
    cursor: &mut Cursor<'_, S>,
    value: Value,
) -> Result<IndexedData, Diagnostic> {
    let fields = expect_map(value)?;
    let mut encoding = None;
    let mut payload = None;
    for _ in 0..fields {
        let key = read_key(cursor)?;
        let value_start = cursor.offset();
        let value = marker::value(cursor)?;
        match key.as_str() {
            "encoding" => encoding = Some(value_range(cursor, value_start, value)?),
            "data" => payload = Some(binary_range(cursor, value)?),
            _ => marker::skip(cursor, value)?,
        }
    }
    match (encoding, payload) {
        (Some(encoding), Some(payload)) => Ok(IndexedData { encoding, payload }),
        _ => Err(error("BinaryCIF encoded data is incomplete")),
    }
}

fn value_range<S: SourceBytes>(
    cursor: &mut Cursor<'_, S>,
    start: u64,
    value: Value,
) -> Result<Range<u64>, Diagnostic> {
    marker::skip(cursor, value)?;
    Ok(start..cursor.offset())
}

fn binary_range<S: SourceBytes>(
    cursor: &mut Cursor<'_, S>,
    value: Value,
) -> Result<Range<u64>, Diagnostic> {
    let Value::Binary(length) = value else {
        return Err(error("BinaryCIF data payload is not MessagePack binary"));
    };
    let start = cursor.offset();
    cursor.advance(length)?;
    Ok(start..cursor.offset())
}

fn read_key<S: SourceBytes>(cursor: &mut Cursor<'_, S>) -> Result<String, Diagnostic> {
    let value = marker::value(cursor)?;
    read_text(cursor, value)
}

fn read_text<S: SourceBytes>(
    cursor: &mut Cursor<'_, S>,
    value: Value,
) -> Result<String, Diagnostic> {
    let Value::Text(length) = value else {
        return Err(error("BinaryCIF map key or name is not text"));
    };
    cursor.text(length)
}

fn expect_map(value: Value) -> Result<u64, Diagnostic> {
    match value {
        Value::Map(length) => Ok(length),
        _ => Err(error("BinaryCIF container value is not a map")),
    }
}

fn expect_array(value: Value) -> Result<u64, Diagnostic> {
    match value {
        Value::Array(length) => Ok(length),
        _ => Err(error("BinaryCIF container value is not an array")),
    }
}

fn expect_unsigned(value: Value) -> Result<u64, Diagnostic> {
    match value {
        Value::Unsigned(number) => Ok(number),
        _ => Err(error("BinaryCIF row count is not an unsigned integer")),
    }
}

fn index_memory(error: impl std::fmt::Display) -> Diagnostic {
    Diagnostic::new(pdbiox_core::Code::E1902).with_context("reason", error.to_string())
}

#[cfg(test)]
#[path = "index_tests.rs"]
mod tests;
