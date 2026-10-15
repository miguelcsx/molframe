//! Value and error helpers shared by the owning and borrowing container readers.

use crate::codec::{Decoded, length_error};
use molframe_cif::CifValue;
use molframe_core::diagnostic::{Code, Diagnostic};
use std::sync::Arc;

pub(crate) fn masked_value(
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

pub(crate) fn container_error(error: &rmp_serde::decode::Error) -> Diagnostic {
    let message = error.to_string();
    let code = if message.contains("unknown variant") {
        Code::E1402
    } else {
        Code::E1401
    };
    Diagnostic::new(code).with_message(message)
}
