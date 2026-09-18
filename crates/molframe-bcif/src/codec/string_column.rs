//! Compact decoded string columns with shared dictionary values.

use indexmap::IndexMap;
use molframe_core::diagnostic::{Code, Diagnostic};
use std::iter::FusedIterator;
use std::ops::Range;
use std::sync::Arc;

/// A decoded string column stored as unique shared values and one `u32` per row.
///
/// Storage grows with the unique UTF-8 data plus four bytes per row. Cloning the
/// column shares both arrays, and cloning a value for an owning DOM only bumps
/// the corresponding [`Arc`] reference count.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DecodedStringColumn {
    dictionary: Arc<[Arc<str>]>,
    indices: Arc<[u32]>,
}

impl DecodedStringColumn {
    /// Builds a column and canonicalises duplicate dictionary values.
    ///
    /// `indices` address the supplied dictionary before canonicalisation.
    ///
    /// # Errors
    ///
    /// Returns `E1401` when an index is outside the supplied dictionary or the
    /// number of unique values exceeds the `u32` representation.
    pub fn new(dictionary: Vec<Arc<str>>, indices: Vec<u32>) -> Result<Self, Diagnostic> {
        let mut lookup = IndexMap::<Arc<str>, u32>::with_capacity(dictionary.len());
        let mut unique = Vec::with_capacity(dictionary.len());
        let mut remap = Vec::with_capacity(dictionary.len());
        for value in dictionary {
            let compact = if let Some(index) = lookup.get(value.as_ref()) {
                *index
            } else {
                let Ok(index) = u32::try_from(unique.len()) else {
                    return Err(dictionary_too_large());
                };
                lookup.insert(Arc::clone(&value), index);
                unique.push(value);
                index
            };
            remap.push(compact);
        }

        let mut compact_indices = Vec::with_capacity(indices.len());
        for (row, index) in indices.into_iter().enumerate() {
            let Ok(source) = usize::try_from(index) else {
                return Err(invalid_index(row, index, remap.len()));
            };
            let Some(compact) = remap.get(source) else {
                return Err(invalid_index(row, index, remap.len()));
            };
            compact_indices.push(*compact);
        }
        Ok(Self::from_validated_parts(unique, compact_indices))
    }

    pub(crate) fn from_validated_parts(dictionary: Vec<Arc<str>>, indices: Vec<u32>) -> Self {
        Self {
            dictionary: dictionary.into(),
            indices: indices.into(),
        }
    }

    pub(crate) fn into_parts(self) -> (Arc<[Arc<str>]>, Arc<[u32]>) {
        (self.dictionary, self.indices)
    }

    /// Number of rows.
    #[must_use]
    pub fn len(&self) -> usize {
        self.indices.len()
    }

    /// Whether the column has no rows.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.indices.is_empty()
    }

    /// Unique shared values in first-seen order.
    #[must_use]
    pub fn dictionary(&self) -> &[Arc<str>] {
        &self.dictionary
    }

    /// One dictionary index per row.
    #[must_use]
    pub fn indices(&self) -> &[u32] {
        &self.indices
    }

    /// A row value as borrowed text.
    #[must_use]
    pub fn get(&self, row: usize) -> Option<&str> {
        self.get_shared(row).map(AsRef::as_ref)
    }

    /// A row value in its shared owning representation.
    #[must_use]
    pub fn get_shared(&self, row: usize) -> Option<&Arc<str>> {
        let index = *self.indices.get(row)?;
        let index = usize::try_from(index).ok()?;
        self.dictionary.get(index)
    }

    /// Iterates over borrowed row values without materialising strings.
    #[must_use]
    pub fn iter(&self) -> DecodedStringIter<'_> {
        DecodedStringIter {
            column: self,
            rows: 0..self.len(),
        }
    }
}

/// Borrowing iterator over a compact decoded string column.
#[derive(Clone, Debug)]
pub struct DecodedStringIter<'a> {
    column: &'a DecodedStringColumn,
    rows: Range<usize>,
}

impl<'a> Iterator for DecodedStringIter<'a> {
    type Item = &'a str;

    fn next(&mut self) -> Option<Self::Item> {
        self.column.get(self.rows.next()?)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.rows.size_hint()
    }
}

impl DoubleEndedIterator for DecodedStringIter<'_> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.column.get(self.rows.next_back()?)
    }
}

impl ExactSizeIterator for DecodedStringIter<'_> {
    fn len(&self) -> usize {
        self.rows.len()
    }
}

impl FusedIterator for DecodedStringIter<'_> {}

impl<'a> IntoIterator for &'a DecodedStringColumn {
    type Item = &'a str;
    type IntoIter = DecodedStringIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

fn invalid_index(row: usize, index: u32, dictionary_len: usize) -> Diagnostic {
    Diagnostic::new(Code::E1401)
        .with_message("decoded string index is outside the dictionary")
        .with_context("row", row.to_string())
        .with_context("index", index.to_string())
        .with_context("dictionary length", dictionary_len.to_string())
}

fn dictionary_too_large() -> Diagnostic {
    Diagnostic::new(Code::E1401).with_message("decoded string dictionary exceeds u32 indices")
}
