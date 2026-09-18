//! `MessagePack` marker decoding and allocation-free value skipping.

use super::cursor::{Cursor, error};
use molframe_core::{Diagnostic, SourceBytes};

#[derive(Clone, Copy)]
pub(super) enum Value {
    Map(u64),
    Array(u64),
    Text(u64),
    Binary(u64),
    Unsigned(u64),
    Nil,
    Other,
}

pub(super) fn value<S: SourceBytes>(cursor: &mut Cursor<'_, S>) -> Result<Value, Diagnostic> {
    let marker = cursor.byte()?;
    match marker {
        0x00..=0x7f => Ok(Value::Unsigned(u64::from(marker))),
        0x80..=0x8f => Ok(Value::Map(u64::from(marker & 0x0f))),
        0x90..=0x9f => Ok(Value::Array(u64::from(marker & 0x0f))),
        0xa0..=0xbf => Ok(Value::Text(u64::from(marker & 0x1f))),
        0xc0 => Ok(Value::Nil),
        0xc2 | 0xc3 | 0xe0..=0xff => Ok(Value::Other),
        0xc4 => binary(cursor, 1),
        0xc5 => binary(cursor, 2),
        0xc6 => binary(cursor, 4),
        0xc7 => extension(cursor, 1),
        0xc8 => extension(cursor, 2),
        0xc9 => extension(cursor, 4),
        0xca | 0xd2 => fixed(cursor, 4),
        0xcb | 0xd3 => fixed(cursor, 8),
        0xcc => unsigned(cursor, 1),
        0xcd => unsigned(cursor, 2),
        0xce => unsigned(cursor, 4),
        0xcf => unsigned(cursor, 8),
        0xd0 => fixed(cursor, 1),
        0xd1 => fixed(cursor, 2),
        0xd4 => extension_fixed(cursor, 1),
        0xd5 => extension_fixed(cursor, 2),
        0xd6 => extension_fixed(cursor, 4),
        0xd7 => extension_fixed(cursor, 8),
        0xd8 => extension_fixed(cursor, 16),
        0xd9 => text(cursor, 1),
        0xda => text(cursor, 2),
        0xdb => text(cursor, 4),
        0xdc => array(cursor, 2),
        0xdd => array(cursor, 4),
        0xde => map(cursor, 2),
        0xdf => map(cursor, 4),
        _ => Err(error(format!("reserved MessagePack marker 0x{marker:02x}"))),
    }
}

pub(super) fn skip<S: SourceBytes>(
    cursor: &mut Cursor<'_, S>,
    first: Value,
) -> Result<(), Diagnostic> {
    let mut pending = 1_u64;
    let mut current = Some(first);
    while pending > 0 {
        let item = match current.take() {
            Some(item) => item,
            None => value(cursor)?,
        };
        pending -= 1;
        match item {
            Value::Map(length) => pending = add_pending(pending, length.saturating_mul(2))?,
            Value::Array(length) => pending = add_pending(pending, length)?,
            Value::Text(length) | Value::Binary(length) => cursor.advance(length)?,
            Value::Unsigned(_) | Value::Nil | Value::Other => {}
        }
    }
    Ok(())
}

fn binary<S: SourceBytes>(cursor: &mut Cursor<'_, S>, bytes: usize) -> Result<Value, Diagnostic> {
    Ok(Value::Binary(cursor.unsigned(bytes)?))
}

fn text<S: SourceBytes>(cursor: &mut Cursor<'_, S>, bytes: usize) -> Result<Value, Diagnostic> {
    Ok(Value::Text(cursor.unsigned(bytes)?))
}

fn unsigned<S: SourceBytes>(cursor: &mut Cursor<'_, S>, bytes: usize) -> Result<Value, Diagnostic> {
    Ok(Value::Unsigned(cursor.unsigned(bytes)?))
}

fn array<S: SourceBytes>(cursor: &mut Cursor<'_, S>, bytes: usize) -> Result<Value, Diagnostic> {
    Ok(Value::Array(cursor.unsigned(bytes)?))
}

fn map<S: SourceBytes>(cursor: &mut Cursor<'_, S>, bytes: usize) -> Result<Value, Diagnostic> {
    Ok(Value::Map(cursor.unsigned(bytes)?))
}

fn fixed<S: SourceBytes>(cursor: &mut Cursor<'_, S>, bytes: u64) -> Result<Value, Diagnostic> {
    cursor.advance(bytes)?;
    Ok(Value::Other)
}

fn extension<S: SourceBytes>(
    cursor: &mut Cursor<'_, S>,
    bytes: usize,
) -> Result<Value, Diagnostic> {
    let length = cursor.unsigned(bytes)?;
    cursor.advance(length.saturating_add(1))?;
    Ok(Value::Other)
}

fn extension_fixed<S: SourceBytes>(
    cursor: &mut Cursor<'_, S>,
    bytes: u64,
) -> Result<Value, Diagnostic> {
    cursor.advance(bytes.saturating_add(1))?;
    Ok(Value::Other)
}

fn add_pending(pending: u64, extra: u64) -> Result<u64, Diagnostic> {
    pending
        .checked_add(extra)
        .ok_or_else(|| error("MessagePack nesting count overflow"))
}
