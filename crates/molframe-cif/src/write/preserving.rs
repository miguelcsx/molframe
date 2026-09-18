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
use molframe_core::io::TextOutput;
use std::fmt::{self, Display, Formatter, Write as _};
use std::io::{self, Write as IoWrite};

/// Writes a document back out, preserving what it held.
///
/// Categories, items and their order survive, and so does the quoting of every
/// value, so a category this library has no interpretation for comes back
/// unharmed.
#[must_use]
pub fn write_preserving(document: &Document) -> String {
    let mut out = String::new();
    render_preserving(document, &mut out);
    out
}

/// Streams a document while preserving its retained category/value order.
///
/// # Errors
///
/// Returns the first destination I/O error.
pub fn write_preserving_to<W: IoWrite>(document: &Document, output: &mut W) -> io::Result<()> {
    let mut output = TextOutput::new(output);
    render_preserving(document, &mut output);
    output.finish()
}

fn render_preserving(document: &Document, out: &mut impl fmt::Write) {
    for block in document.blocks() {
        let _ = writeln!(out, "data_{}", block.name());
        for category in block.categories() {
            write_category(out, category);
        }
    }
}

fn write_category(out: &mut impl fmt::Write, category: &crate::document::Category) {
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
                RenderedValue { value, quoting }
            );
        }
        let _ = out.write_str("#\n");
        return;
    }

    let _ = out.write_str("loop_\n");
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
                let _ = out.write_char(' ');
            }
            first = false;
            let _ = write!(out, "{}", RenderedValue { value, quoting });
        }
        let _ = out.write_char('\n');
    }
    let _ = out.write_str("#\n");
}

struct RenderedValue<'a> {
    value: &'a CifValue,
    quoting: Option<Quoting>,
}

impl Display for RenderedValue<'_> {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        let text = match self.value {
            CifValue::Inapplicable => return formatter.write_char('.'),
            CifValue::Unknown => return formatter.write_char('?'),
            CifValue::Integer(number) => return Display::fmt(number, formatter),
            CifValue::Float(number) => return Display::fmt(number, formatter),
            CifValue::Text(text) => text,
        };
        match self.quoting {
            Some(Quoting::Text) => write!(formatter, "\n;{text}\n;"),
            Some(Quoting::Double) => write!(formatter, "\"{text}\""),
            Some(Quoting::Single) => write!(formatter, "'{text}'"),
            _ => Display::fmt(&super::value::quoted(text), formatter),
        }
    }
}

#[cfg(test)]
#[path = "preserving_tests.rs"]
mod tests;
