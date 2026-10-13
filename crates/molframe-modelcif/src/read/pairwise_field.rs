//! Owned, reusable copy of one CIF scalar for the streaming PAE reader.
//!
//! The streaming reader keeps at most one pending table row, so each field
//! reuses its text buffer instead of allocating per cell.

use molframe_cif::CifScalar;

/// Largest integer a binary64 value holds exactly.
const EXACT_INTEGER: i64 = 1 << 53;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
enum Kind {
    #[default]
    Missing,
    Inapplicable,
    Unknown,
    Text,
    Integer,
    Float,
}

#[derive(Debug, Default)]
pub(super) struct Field {
    kind: Kind,
    text: String,
    integer: i64,
    float: f64,
}

impl Field {
    pub(super) fn set(&mut self, value: CifScalar<'_>) {
        match value {
            CifScalar::Inapplicable => self.kind = Kind::Inapplicable,
            CifScalar::Unknown => self.kind = Kind::Unknown,
            CifScalar::Text(text) => {
                self.text.clear();
                self.text.push_str(text);
                self.kind = Kind::Text;
            }
            CifScalar::Integer(value) => {
                self.integer = value;
                self.kind = Kind::Integer;
            }
            CifScalar::Float(value) => {
                self.float = value;
                self.kind = Kind::Float;
            }
        }
    }

    pub(super) fn copy_from(&mut self, other: &Self) {
        self.kind = other.kind;
        self.integer = other.integer;
        self.float = other.float;
        self.text.clear();
        self.text.push_str(&other.text);
    }

    pub(super) fn scalar(&self) -> CifScalar<'_> {
        match self.kind {
            Kind::Missing | Kind::Unknown => CifScalar::Unknown,
            Kind::Inapplicable => CifScalar::Inapplicable,
            Kind::Text => CifScalar::Text(&self.text),
            Kind::Integer => CifScalar::Integer(self.integer),
            Kind::Float => CifScalar::Float(self.float),
        }
    }

    pub(super) fn same(&self, other: &Self) -> bool {
        self.kind == other.kind
            && match self.kind {
                Kind::Text => self.text == other.text,
                Kind::Integer => self.integer == other.integer,
                Kind::Float => self.float.to_bits() == other.float.to_bits(),
                _ => true,
            }
    }

    pub(super) fn text(&self) -> Option<&str> {
        (self.kind == Kind::Text).then_some(self.text.as_str())
    }

    /// Text of an identifier written as a word or as a whole number.
    pub(super) fn identifier(&self) -> Option<String> {
        match self.kind {
            Kind::Text => Some(self.text.clone()),
            Kind::Integer => Some(self.integer.to_string()),
            _ => None,
        }
    }

    pub(super) fn integer(&self) -> Option<i64> {
        (self.kind == Kind::Integer).then_some(self.integer)
    }

    pub(super) fn is_float(&self) -> bool {
        self.kind == Kind::Float
    }

    /// Finite numeric value, when binary64 holds it exactly.
    pub(super) fn number(&self) -> Option<f64> {
        match self.kind {
            #[allow(clippy::cast_precision_loss)]
            Kind::Integer if self.integer.abs() <= EXACT_INTEGER => Some(self.integer as f64),
            Kind::Float if self.float.is_finite() => Some(self.float),
            _ => None,
        }
    }
}
