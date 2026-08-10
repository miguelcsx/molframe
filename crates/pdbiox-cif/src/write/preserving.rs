//! Writing, under two contracts that are stated rather than implied.
//!
//! **Canonical** produces a valid, self-consistent file. It does not promise to
//! reproduce the original's formatting, and it does not carry categories the
//! structure does not model.
//!
//! **Preserving** writes a document back out with its categories, its items and
//! their order intact, and with every untouched value written the way it was
//! found — including whether it was quoted.
//!
//! Neither promises a byte-exact round trip after the data has been changed.
//! That promise cannot be kept honestly once a value differs in width, and a
//! library that implies it will be trusted where it should not be.

use crate::document::{CifValue, Document};
use crate::lexer::Quoting;
use std::fmt::Write as _;

/// Writes a document back out, preserving what it held.
///
/// Categories, items and their order survive, and so does the quoting of every
/// value, so a category this library has no interpretation for comes back
/// unharmed.
#[must_use]
pub fn write_preserving(document: &Document) -> String {
    let mut out = String::new();
    for block in document.blocks() {
        let _ = writeln!(out, "data_{}", block.name());
        for category in block.categories() {
            write_category(&mut out, category);
        }
    }
    out
}

fn write_category(out: &mut String, category: &crate::document::Category) {
    let rows = category.row_count();
    if rows == 1 {
        for item in category.items() {
            let Some(value) = category.value(item, 0) else {
                continue;
            };
            let quoting = category.column(item).and_then(|column| column.quoting(0));
            let _ = writeln!(
                out,
                "_{}.{:<30} {}",
                category.name(),
                item,
                render(value, quoting)
            );
        }
        out.push_str("#\n");
        return;
    }

    out.push_str("loop_\n");
    for item in category.items() {
        let _ = writeln!(out, "_{}.{item}", category.name());
    }
    for row in 0..rows {
        let mut first = true;
        for item in category.items() {
            let Some(value) = category.value(item, row) else {
                continue;
            };
            let quoting = category.column(item).and_then(|column| column.quoting(row));
            if !first {
                out.push(' ');
            }
            first = false;
            out.push_str(&render(value, quoting));
        }
        out.push('\n');
    }
    out.push_str("#\n");
}

/// Renders one value the way it was written, quoting it if it needs quoting.
fn render(value: &CifValue, quoting: Option<Quoting>) -> String {
    let text = match value {
        CifValue::Inapplicable => return ".".to_owned(),
        CifValue::Unknown => return "?".to_owned(),
        CifValue::Integer(number) => return number.to_string(),
        CifValue::Float(number) => return format!("{number}"),
        CifValue::Text(text) => text,
    };
    match quoting {
        Some(Quoting::Text) => format!("\n;{text}\n;"),
        Some(Quoting::Double) => format!("\"{text}\""),
        Some(Quoting::Single) => format!("'{text}'"),
        _ => super::value::quote_text(text),
    }
}

#[cfg(test)]
#[path = "preserving_tests.rs"]
mod tests;
