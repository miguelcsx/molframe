//! Fixed-width borrowed `atom_site` rows.

use crate::document::CifValueRef;
use crate::lower::{AtomSiteRow, Field};
use std::borrow::Cow;

const FIELD_COUNT: usize = Field::COUNT;

pub(super) struct AtomRow<'input> {
    values: [Option<CifValueRef<'input>>; FIELD_COUNT],
    row: usize,
}

impl<'input> AtomRow<'input> {
    pub(super) const fn new() -> Self {
        Self {
            values: [None; FIELD_COUNT],
            row: 0,
        }
    }

    pub(super) fn set(&mut self, field: Field, value: CifValueRef<'input>) {
        if let Some(slot) = self.values.get_mut(field.position()) {
            *slot = Some(value);
        }
    }

    pub(super) fn supports(item: &str) -> bool {
        Field::from_item(item).is_some()
    }

    /// The value at `field`, by position rather than by name comparison.
    fn value(&self, field: Field) -> Option<CifValueRef<'input>> {
        self.values.get(field.position()).copied().flatten()
    }

    pub(super) fn advance(&mut self) {
        self.values.fill(None);
        self.row += 1;
    }
}

impl AtomSiteRow for AtomRow<'_> {
    fn row(&self) -> usize {
        self.row
    }

    fn text(&self, field: Field) -> Option<&str> {
        self.value(field)?.as_str()
    }

    fn identifier(&self, field: Field) -> Option<Cow<'_, str>> {
        self.value(field)?.as_identifier()
    }

    fn integer(&self, field: Field) -> Option<i64> {
        self.value(field)?.as_integer()
    }

    fn float(&self, field: Field) -> Option<f64> {
        self.value(field)?.as_float()
    }

    fn is_recorded(&self, field: Field) -> bool {
        self.value(field).is_some_and(CifValueRef::is_recorded)
    }
}
