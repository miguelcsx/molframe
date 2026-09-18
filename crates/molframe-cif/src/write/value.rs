//! One canonical renderer for CIF values and identifiers.

use crate::document::CifValue;
use std::fmt::{self, Display, Formatter, Write};

/// Renders a typed value as one CIF token or text field.
#[must_use]
pub fn render_value(value: &CifValue) -> String {
    match value {
        CifValue::Inapplicable => ".".to_owned(),
        CifValue::Unknown => "?".to_owned(),
        CifValue::Integer(value) => value.to_string(),
        CifValue::Float(value) => value.to_string(),
        CifValue::Text(value) => quote_text(value),
    }
}

/// Quotes arbitrary text without changing its value when reparsed.
#[must_use]
pub fn quote_text(text: &str) -> String {
    let mut rendered = String::new();
    write_quoted(&mut rendered, text);
    rendered
}

/// Appends arbitrary text as one CIF token or text field without a temporary allocation.
///
/// This is the streaming counterpart of [`quote_text`]. It uses exactly the
/// same canonical quoting decisions while allowing large columnar writers to
/// reuse their destination buffer.
pub fn write_quoted(output: &mut impl Write, text: &str) {
    let _ = write_quoted_to(output, text);
}

pub(crate) struct Quoted<'a>(&'a str);

pub(crate) const fn quoted(text: &str) -> Quoted<'_> {
    Quoted(text)
}

impl Display for Quoted<'_> {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        write_quoted_to(formatter, self.0)
    }
}

fn write_quoted_to(output: &mut impl Write, text: &str) -> fmt::Result {
    if text.contains('\n') || (text.contains('\'') && text.contains('"')) {
        output.write_str("\n;")?;
        output.write_str(text)?;
        return output.write_str("\n;");
    }
    if !needs_quoting(text) {
        return output.write_str(text);
    }
    if text.contains('\'') {
        output.write_char('"')?;
        output.write_str(text)?;
        output.write_char('"')
    } else {
        output.write_char('\'')?;
        output.write_str(text)?;
        output.write_char('\'')
    }
}

pub(crate) fn needs_quoting(text: &str) -> bool {
    text.is_empty()
        || text.chars().any(char::is_whitespace)
        || matches!(text, "." | "?")
        || text.starts_with(['_', '#', '\'', '"', '[', ']', '$', ';'])
        || ["loop_", "stop_", "global_"]
            .iter()
            .any(|keyword| text.eq_ignore_ascii_case(keyword))
        || text.get(..5).is_some_and(|prefix| {
            prefix.eq_ignore_ascii_case("data_") || prefix.eq_ignore_ascii_case("save_")
        })
}

#[cfg(test)]
#[path = "value_tests.rs"]
mod tests;
