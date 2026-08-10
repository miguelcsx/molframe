//! One canonical renderer for CIF values and identifiers.

use crate::document::CifValue;

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
    if text.contains('\n') || (text.contains('\'') && text.contains('"')) {
        return format!("\n;{text}\n;");
    }
    if !needs_quoting(text) {
        return text.to_owned();
    }
    if text.contains('\'') {
        format!("\"{text}\"")
    } else {
        format!("'{text}'")
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
