//! The lossless layer: what the file said, before anything interprets it.
//!
//! A document keeps every category, every item and the order they appeared in,
//! including the ones pdbiox has no interpretation for. That is what lets a file
//! survive a read and a write even where this library predates the dictionary
//! extension it uses, and it is the only layer that can honestly claim to
//! preserve what was written.
//!
//! Values keep the distinction the format draws between "does not apply", "not
//! recorded" and "here is the value", because those mean different things and
//! collapsing them into one absent value throws away a statement the depositor
//! made deliberately.

use crate::lexer::Quoting;
use indexmap::IndexMap;
use num_traits::ToPrimitive;
use pdbiox_core::span::ByteSpan;
use std::borrow::Cow;

/// One value of one item.
#[derive(Clone, PartialEq, Debug)]
pub enum CifValue {
    /// Written `.` — the item does not apply here.
    Inapplicable,
    /// Written `?` — the value exists but was not recorded.
    Unknown,
    /// Text.
    Text(Box<str>),
    /// A whole number.
    ///
    /// Kept as an integer rather than as a float: a serial number past two to
    /// the fifty-third would otherwise stop round-tripping exactly, and the
    /// files that reach those numbers are precisely the large ones.
    Integer(i64),
    /// A number with a fractional part.
    Float(f64),
}

impl CifValue {
    /// Interprets a lexed value, deciding what kind of thing it is.
    ///
    /// Only a bare value can be a number or a sentinel. `'.'` written in quotes
    /// is the one-character string, not the sentinel, and the format is explicit
    /// about the difference.
    #[must_use]
    pub fn parse(text: &str, quoting: Quoting) -> Self {
        if quoting != Quoting::Bare {
            return Self::Text(text.into());
        }
        match text {
            "." => Self::Inapplicable,
            "?" => Self::Unknown,
            _ => {
                if let Ok(value) = text.parse::<i64>() {
                    return Self::Integer(value);
                }
                if let Ok(value) = text.parse::<f64>() {
                    return Self::Float(value);
                }
                Self::Text(text.into())
            }
        }
    }

    /// The value as text, or `None` when it is a sentinel.
    #[must_use]
    pub fn as_str(&self) -> Option<&str> {
        match self {
            Self::Text(text) => Some(text),
            _ => None,
        }
    }

    /// The value as an identifier, whatever kind of thing it was written as.
    ///
    /// Identifiers in this format are frequently numeric — entity numbers and
    /// chain labels especially — and reading them only as text drops them
    /// silently. Nothing is allocated unless the value really was a number.
    #[must_use]
    pub fn as_identifier(&self) -> Option<Cow<'_, str>> {
        match self {
            Self::Text(text) => Some(Cow::Borrowed(text)),
            Self::Integer(value) => Some(Cow::Owned(value.to_string())),
            Self::Float(value) => Some(Cow::Owned(value.to_string())),
            Self::Inapplicable | Self::Unknown => None,
        }
    }

    /// The value as a whole number, where it is one.
    #[must_use]
    pub fn as_integer(&self) -> Option<i64> {
        match self {
            Self::Integer(value) => Some(*value),
            _ => None,
        }
    }

    /// The value as a number, whether it was written with a fractional part.
    #[must_use]
    pub fn as_float(&self) -> Option<f64> {
        match self {
            Self::Integer(value) => value.to_f64(),
            Self::Float(value) => Some(*value),
            _ => None,
        }
    }

    /// Returns true for a value that was recorded.
    #[must_use]
    pub const fn is_recorded(&self) -> bool {
        !matches!(self, Self::Inapplicable | Self::Unknown)
    }
}

/// One item's values, in row order.
#[derive(Clone, Debug, Default)]
pub struct Column {
    values: Vec<CifValue>,
    quoting: Vec<Quoting>,
}

impl Column {
    /// Appends a value and how it was written.
    pub fn push(&mut self, value: CifValue, quoting: Quoting) {
        self.values.push(value);
        self.quoting.push(quoting);
    }

    /// The number of rows.
    #[must_use]
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// Returns true when the column holds nothing.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// The value in one row.
    #[must_use]
    pub fn get(&self, row: usize) -> Option<&CifValue> {
        self.values.get(row)
    }

    /// How the value in one row was written, so a preserving write can put it
    /// back the way it found it.
    #[must_use]
    pub fn quoting(&self, row: usize) -> Option<Quoting> {
        self.quoting.get(row).copied()
    }

    /// Every value, in row order.
    pub fn iter(&self) -> impl Iterator<Item = &CifValue> {
        self.values.iter()
    }
}

/// One category and its items.
#[derive(Clone, Debug)]
pub struct Category {
    name: Box<str>,
    columns: IndexMap<Box<str>, Column>,
    span: ByteSpan,
}

impl Category {
    /// Creates an empty category.
    #[must_use]
    pub fn new(name: impl Into<Box<str>>, span: ByteSpan) -> Self {
        Self {
            name: name.into(),
            columns: IndexMap::new(),
            span,
        }
    }

    /// The category's name, without a leading underscore.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Where the category was found.
    #[must_use]
    pub const fn span(&self) -> ByteSpan {
        self.span
    }

    /// The number of rows, taken from the longest item.
    #[must_use]
    pub fn row_count(&self) -> usize {
        // A category with no items has no rows, rather than an unknown number.
        match self.columns.values().map(Column::len).max() {
            Some(rows) => rows,
            None => 0,
        }
    }

    /// The number of items.
    #[must_use]
    pub fn len(&self) -> usize {
        self.columns.len()
    }

    /// Returns true when the category holds no items.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.columns.is_empty()
    }

    /// One item's values.
    #[must_use]
    pub fn column(&self, item: &str) -> Option<&Column> {
        self.columns.get(item)
    }

    /// One item's values, creating the item if it is new.
    pub fn column_mut(&mut self, item: &str) -> &mut Column {
        self.columns.entry(item.into()).or_default()
    }

    /// The item names, in the order the file gave them.
    pub fn items(&self) -> impl Iterator<Item = &str> {
        self.columns.keys().map(Box::as_ref)
    }

    /// One value, by item and row.
    #[must_use]
    pub fn value(&self, item: &str, row: usize) -> Option<&CifValue> {
        self.column(item)?.get(row)
    }

    /// One value as text, by item and row.
    #[must_use]
    pub fn text(&self, item: &str, row: usize) -> Option<&str> {
        self.value(item, row)?.as_str()
    }

    /// One value as an identifier, by item and row.
    #[must_use]
    pub fn identifier(&self, item: &str, row: usize) -> Option<Cow<'_, str>> {
        self.value(item, row)?.as_identifier()
    }
}

/// One data block.
#[derive(Clone, Debug)]
pub struct DataBlock {
    name: Box<str>,
    categories: IndexMap<Box<str>, Category>,
}

impl DataBlock {
    /// Creates an empty block.
    #[must_use]
    pub fn new(name: impl Into<Box<str>>) -> Self {
        Self {
            name: name.into(),
            categories: IndexMap::new(),
        }
    }

    /// The block's name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// One category.
    #[must_use]
    pub fn category(&self, name: &str) -> Option<&Category> {
        self.categories.get(name)
    }

    /// One category, creating it if it is new.
    pub fn category_mut(&mut self, name: &str, span: ByteSpan) -> &mut Category {
        self.categories
            .entry(name.into())
            .or_insert_with(|| Category::new(name, span))
    }

    /// The categories, in the order the file gave them.
    pub fn categories(&self) -> impl Iterator<Item = &Category> {
        self.categories.values()
    }

    /// The number of categories.
    #[must_use]
    pub fn len(&self) -> usize {
        self.categories.len()
    }

    /// Returns true when the block holds no categories.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.categories.is_empty()
    }
}

/// A whole file, losslessly.
///
/// # Examples
///
/// ```
/// use pdbiox_cif::{Document, parse};
/// use pdbiox_core::io::InputBuffer;
///
/// let input = InputBuffer::from_bytes(b"data_test\n_entry.id 1ABC\n".to_vec());
/// let (document, _) = parse(&input)?;
///
/// let id = document.first_block().and_then(|block| block.category("entry"))
///     .and_then(|entry| entry.text("id", 0));
/// assert_eq!(id, Some("1ABC"));
/// # Ok::<(), Vec<pdbiox_core::diagnostic::Diagnostic>>(())
/// ```
#[derive(Clone, Debug, Default)]
pub struct Document {
    blocks: Vec<DataBlock>,
}

impl Document {
    /// Creates an empty document.
    #[must_use]
    pub const fn new() -> Self {
        Self { blocks: Vec::new() }
    }

    /// Appends a block.
    pub fn push(&mut self, block: DataBlock) {
        self.blocks.push(block);
    }

    /// The first block, which is the one a single-entry file has.
    #[must_use]
    pub fn first_block(&self) -> Option<&DataBlock> {
        self.blocks.first()
    }

    /// The last block, which is the one still being filled while parsing.
    pub fn last_block_mut(&mut self) -> Option<&mut DataBlock> {
        self.blocks.last_mut()
    }

    /// Every block, in order.
    pub fn blocks(&self) -> impl Iterator<Item = &DataBlock> {
        self.blocks.iter()
    }

    /// The number of blocks.
    #[must_use]
    pub fn len(&self) -> usize {
        self.blocks.len()
    }

    /// Returns true when the document holds no blocks.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.blocks.is_empty()
    }
}

#[cfg(test)]
#[path = "model_tests.rs"]
mod tests;
