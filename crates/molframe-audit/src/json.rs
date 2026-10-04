//! A deterministic JSON value, written in the order it was built.

use molframe_core::contract::push_json_string;
use std::fmt::Write;

/// A JSON value whose objects keep the order their keys were added in.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Json {
    Null,
    Bool(bool),
    Number(f64),
    Text(String),
    Array(Vec<Self>),
    Object(Vec<(&'static str, Self)>),
}

impl Json {
    pub(crate) fn text(value: impl Into<String>) -> Self {
        Self::Text(value.into())
    }

    /// A reference to another entity of the graph.
    pub(crate) fn reference(id: &str) -> Self {
        Self::Object(vec![("@id", Self::text(id))])
    }

    pub(crate) fn render(&self) -> String {
        let mut output = String::new();
        self.write(&mut output);
        output
    }

    fn write(&self, output: &mut String) {
        match self {
            Self::Bool(value) => output.push_str(if *value { "true" } else { "false" }),
            Self::Number(value) if value.is_finite() => {
                let _ = write!(output, "{value}");
            }
            // JSON has no spelling for a non-finite number; the absence is the honest one.
            Self::Null | Self::Number(_) => output.push_str("null"),
            Self::Text(value) => push_json_string(output, value),
            Self::Array(items) => {
                output.push('[');
                for (position, item) in items.iter().enumerate() {
                    if position != 0 {
                        output.push(',');
                    }
                    item.write(output);
                }
                output.push(']');
            }
            Self::Object(members) => {
                output.push('{');
                for (position, (key, value)) in members.iter().enumerate() {
                    if position != 0 {
                        output.push(',');
                    }
                    push_json_string(output, key);
                    output.push(':');
                    value.write(output);
                }
                output.push('}');
            }
        }
    }
}
