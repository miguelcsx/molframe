//! Compact columnar retention for forward-compatible `ma_*` categories.
//!
//! Lookup is `O(1)` once a typed iterator resolves its columns. Identifier
//! columns dictionary-encode repeated keys, integer columns use the narrowest
//! exact signed representation, and numeric columns retain integer-versus-float
//! semantics without one enum allocation per cell.

use super::packed::{PackedIndices, PackedIntegers, reserved};
use super::schema::identifier_item;
use crate::{MemoryBudget, ModelCifError};
use pdbiox_cif::{Category, CifValue};
use std::collections::HashMap;
use std::mem::size_of;
use std::sync::Arc;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum ValueRef<'a> {
    Inapplicable,
    Unknown,
    Text(&'a Arc<str>),
    Integer(i64),
    Float(f64),
}

impl ValueRef<'_> {
    pub(crate) const fn integer(self) -> Option<i64> {
        match self {
            Self::Integer(value) => Some(value),
            _ => None,
        }
    }

    pub(crate) fn number(self) -> Option<f64> {
        match self {
            Self::Integer(value) => num_traits::ToPrimitive::to_f64(&value),
            Self::Float(value) => Some(value),
            Self::Inapplicable | Self::Unknown | Self::Text(_) => None,
        }
    }

    pub(crate) fn to_owned(self) -> CifValue {
        match self {
            Self::Inapplicable => CifValue::Inapplicable,
            Self::Unknown => CifValue::Unknown,
            Self::Text(value) => CifValue::Text(Arc::clone(value)),
            Self::Integer(value) => CifValue::Integer(value),
            Self::Float(value) => CifValue::Float(value),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct DictionaryEntry {
    pub(crate) value: CifValue,
    pub(crate) numeric_text: Option<Arc<str>>,
}

impl DictionaryEntry {
    pub(crate) fn new(value: &CifValue) -> Self {
        let numeric_text = match value {
            CifValue::Integer(value) => Some(Arc::from(value.to_string())),
            CifValue::Float(value) => Some(Arc::from(value.to_string())),
            CifValue::Inapplicable | CifValue::Unknown | CifValue::Text(_) => None,
        };
        Self {
            value: value.clone(),
            numeric_text,
        }
    }

    fn identifier(&self) -> Option<&str> {
        match &self.value {
            CifValue::Text(value) => Some(value),
            CifValue::Integer(_) | CifValue::Float(_) => self.numeric_text.as_deref(),
            CifValue::Inapplicable | CifValue::Unknown => None,
        }
    }

    fn as_ref(&self) -> ValueRef<'_> {
        match &self.value {
            CifValue::Inapplicable => ValueRef::Inapplicable,
            CifValue::Unknown => ValueRef::Unknown,
            CifValue::Text(value) => ValueRef::Text(value),
            CifValue::Integer(value) => ValueRef::Integer(*value),
            CifValue::Float(value) => ValueRef::Float(*value),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct DictionaryColumn {
    pub(crate) entries: Vec<DictionaryEntry>,
    pub(crate) indices: PackedIndices,
}

impl DictionaryColumn {
    fn value(&self, row: usize) -> Option<ValueRef<'_>> {
        let index = usize::try_from(self.indices.get(row)?).ok()?;
        self.entries.get(index).map(DictionaryEntry::as_ref)
    }

    fn identifier(&self, row: usize) -> Option<&str> {
        let index = usize::try_from(self.indices.get(row)?).ok()?;
        self.entries.get(index)?.identifier()
    }

    pub(crate) fn bytes(&self) -> usize {
        let entry_bytes = self
            .entries
            .iter()
            .map(|entry| {
                size_of::<DictionaryEntry>()
                    + entry.numeric_text.as_ref().map_or(0, |value| value.len())
                    + entry.value.as_str().map_or(0, str::len)
            })
            .sum::<usize>();
        self.indices.bytes().saturating_add(entry_bytes)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct IntegerColumn {
    pub(crate) values: PackedIntegers,
    pub(crate) states: Option<Vec<u8>>,
}

impl IntegerColumn {
    fn value(&self, row: usize) -> Option<ValueRef<'_>> {
        match self
            .states
            .as_ref()
            .and_then(|states| states.get(row))
            .copied()
        {
            Some(1) => Some(ValueRef::Unknown),
            Some(2) => Some(ValueRef::Inapplicable),
            Some(_) | None => self.values.get(row).map(ValueRef::Integer),
        }
    }

    fn bytes(&self) -> usize {
        self.values.bytes() + self.states.as_ref().map_or(0, Vec::len)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct NumberColumn {
    pub(crate) payload: Vec<u64>,
    pub(crate) kinds: Option<Vec<u8>>,
}

impl NumberColumn {
    fn value(&self, row: usize) -> Option<ValueRef<'_>> {
        let bits = *self.payload.get(row)?;
        match self
            .kinds
            .as_ref()
            .and_then(|kinds| kinds.get(row))
            .copied()
        {
            Some(1) => Some(ValueRef::Integer(i64::from_ne_bytes(bits.to_ne_bytes()))),
            Some(2) => Some(ValueRef::Unknown),
            Some(3) => Some(ValueRef::Inapplicable),
            Some(_) | None => Some(ValueRef::Float(f64::from_bits(bits))),
        }
    }

    pub(crate) fn bytes(&self) -> usize {
        self.payload.len().saturating_mul(size_of::<u64>())
            + self.kinds.as_ref().map_or(0, Vec::len)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) enum CompactColumn {
    Dictionary(DictionaryColumn),
    Integer(IntegerColumn),
    Number(NumberColumn),
    Dense(Vec<CifValue>),
}

impl CompactColumn {
    fn from_category(category: &Category, item: &str) -> Result<Self, ModelCifError> {
        let force_dictionary = identifier_item(category.name(), item);
        let mut has_text = false;
        let mut has_integer = false;
        let mut has_float = false;
        for row in 0..category.row_count() {
            match category.value(item, row) {
                Some(CifValue::Text(_)) => has_text = true,
                Some(CifValue::Integer(_)) => has_integer = true,
                Some(CifValue::Float(_)) => has_float = true,
                Some(CifValue::Inapplicable | CifValue::Unknown) | None => {}
            }
        }
        if force_dictionary || (has_text && !has_integer && !has_float) {
            return Self::dictionary(category, item);
        }
        if !has_text && !has_float {
            return Self::integers(category, item);
        }
        if !has_text {
            return Self::numbers(category, item);
        }
        Self::dense(category, item)
    }

    fn dictionary(category: &Category, item: &str) -> Result<Self, ModelCifError> {
        let rows = category.row_count();
        let mut entries = Vec::new();
        let mut lookup = HashMap::new();
        lookup
            .try_reserve(rows.min(4_096))
            .map_err(|_| ModelCifError::Allocation)?;
        let mut indices = reserved(rows)?;
        for row in 0..rows {
            let value = value_or_unknown(category, item, row);
            let key = ScalarKey::new(value);
            let index = if let Some(index) = lookup.get(&key).copied() {
                index
            } else {
                let index = u32::try_from(entries.len()).map_err(|_| ModelCifError::Capacity)?;
                entries.push(DictionaryEntry::new(value));
                lookup.insert(key, index);
                index
            };
            indices.push(index);
        }
        let maximum =
            u32::try_from(entries.len().saturating_sub(1)).map_err(|_| ModelCifError::Capacity)?;
        Ok(Self::Dictionary(DictionaryColumn {
            entries,
            indices: PackedIndices::from_u32(indices, maximum)?,
        }))
    }

    fn integers(category: &Category, item: &str) -> Result<Self, ModelCifError> {
        let rows = category.row_count();
        let mut values = reserved(rows)?;
        let mut states = reserved(rows)?;
        let mut has_sentinel = false;
        let mut minimum = 0_i64;
        let mut maximum = 0_i64;
        let mut saw_value = false;
        for row in 0..rows {
            let (value, state) = match category.value(item, row) {
                Some(CifValue::Integer(value)) => (*value, 0),
                Some(CifValue::Inapplicable) => (0, 2),
                Some(CifValue::Unknown) | None => (0, 1),
                Some(CifValue::Text(_) | CifValue::Float(_)) => return Self::dense(category, item),
            };
            if state == 0 {
                minimum = if saw_value { minimum.min(value) } else { value };
                maximum = if saw_value { maximum.max(value) } else { value };
                saw_value = true;
            } else {
                has_sentinel = true;
            }
            values.push(value);
            states.push(state);
        }
        Ok(Self::Integer(IntegerColumn {
            values: PackedIntegers::from_i64(&values, minimum, maximum)?,
            states: has_sentinel.then_some(states),
        }))
    }

    fn numbers(category: &Category, item: &str) -> Result<Self, ModelCifError> {
        let rows = category.row_count();
        let mut payload = reserved(rows)?;
        let mut kinds = reserved(rows)?;
        let mut all_float = true;
        for row in 0..rows {
            let (bits, kind) = match category.value(item, row) {
                Some(CifValue::Float(value)) => (value.to_bits(), 0),
                Some(CifValue::Integer(value)) => (u64::from_ne_bytes(value.to_ne_bytes()), 1),
                Some(CifValue::Unknown) | None => (0, 2),
                Some(CifValue::Inapplicable) => (0, 3),
                Some(CifValue::Text(_)) => return Self::dense(category, item),
            };
            all_float &= kind == 0;
            payload.push(bits);
            kinds.push(kind);
        }
        Ok(Self::Number(NumberColumn {
            payload,
            kinds: (!all_float).then_some(kinds),
        }))
    }

    fn dense(category: &Category, item: &str) -> Result<Self, ModelCifError> {
        let mut values = reserved(category.row_count())?;
        for row in 0..category.row_count() {
            values.push(value_or_unknown(category, item, row).clone());
        }
        Ok(Self::Dense(values))
    }

    pub(crate) fn value(&self, row: usize) -> Option<ValueRef<'_>> {
        match self {
            Self::Dictionary(column) => column.value(row),
            Self::Integer(column) => column.value(row),
            Self::Number(column) => column.value(row),
            Self::Dense(values) => values.get(row).map(|value| match value {
                CifValue::Inapplicable => ValueRef::Inapplicable,
                CifValue::Unknown => ValueRef::Unknown,
                CifValue::Text(value) => ValueRef::Text(value),
                CifValue::Integer(value) => ValueRef::Integer(*value),
                CifValue::Float(value) => ValueRef::Float(*value),
            }),
        }
    }

    pub(crate) fn identifier(&self, row: usize) -> Option<&str> {
        match self {
            Self::Dictionary(column) => column.identifier(row),
            Self::Integer(_) | Self::Number(_) | Self::Dense(_) => None,
        }
    }

    pub(crate) fn bytes(&self) -> usize {
        match self {
            Self::Dictionary(column) => column.bytes(),
            Self::Integer(column) => column.bytes(),
            Self::Number(column) => column.bytes(),
            Self::Dense(values) => values.len().saturating_mul(size_of::<CifValue>()),
        }
    }
}

/// One complete compact `ma_*` category.
#[derive(Clone, Debug, PartialEq)]
pub struct ModelCategory {
    pub(crate) name: Box<str>,
    pub(crate) items: Vec<Box<str>>,
    pub(crate) columns: Vec<CompactColumn>,
    pub(crate) rows: usize,
}

impl ModelCategory {
    pub(crate) fn from_category(
        category: &Category,
        budget: &mut MemoryBudget,
    ) -> Result<Self, ModelCifError> {
        let mut items = Vec::new();
        let mut columns = Vec::new();
        items
            .try_reserve_exact(category.len())
            .map_err(|_| ModelCifError::Allocation)?;
        columns
            .try_reserve_exact(category.len())
            .map_err(|_| ModelCifError::Allocation)?;
        let mut bytes = category.name().len();
        for item in category.items() {
            let column = CompactColumn::from_category(category, item)?;
            bytes = bytes
                .checked_add(item.len())
                .and_then(|value| value.checked_add(column.bytes()))
                .ok_or(ModelCifError::Capacity)?;
            items.push(item.into());
            columns.push(column);
        }
        budget.claim(bytes)?;
        Ok(Self {
            name: category.name().into(),
            items,
            columns,
            rows: category.row_count(),
        })
    }

    /// Category name without a leading underscore.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Dictionary items in source order.
    #[must_use]
    pub fn items(&self) -> &[Box<str>] {
        &self.items
    }

    /// Number of source rows.
    #[must_use]
    pub const fn row_count(&self) -> usize {
        self.rows
    }

    /// Rows in source order, materialised only as borrowed cursors.
    #[must_use]
    pub fn rows(&self) -> ModelRows<'_> {
        ModelRows {
            category: self,
            row: 0,
        }
    }

    /// Looks up one exact CIF value by item and row.
    #[must_use]
    pub fn value(&self, item: &str, row: usize) -> Option<CifValue> {
        self.column(item)?.value(row).map(ValueRef::to_owned)
    }

    pub(crate) fn column(&self, item: &str) -> Option<&CompactColumn> {
        let index = self
            .items
            .iter()
            .position(|candidate| candidate.as_ref() == item)?;
        self.columns.get(index)
    }

    pub(crate) fn value_ref(&self, column: usize, row: usize) -> Option<ValueRef<'_>> {
        self.columns.get(column)?.value(row)
    }
}

/// Borrowed cursor for one compact category row.
#[derive(Clone, Copy, Debug)]
pub struct ModelRow<'a> {
    category: &'a ModelCategory,
    row: usize,
}

impl ModelRow<'_> {
    /// One value in item order.
    #[must_use]
    pub fn value(&self, index: usize) -> Option<CifValue> {
        self.category
            .columns
            .get(index)?
            .value(self.row)
            .map(ValueRef::to_owned)
    }

    pub(crate) fn value_ref(&self, index: usize) -> Option<ValueRef<'_>> {
        self.category.value_ref(index, self.row)
    }
}

/// Allocation-free iteration over compact category rows.
#[derive(Clone, Debug)]
pub struct ModelRows<'a> {
    category: &'a ModelCategory,
    row: usize,
}

impl<'a> Iterator for ModelRows<'a> {
    type Item = ModelRow<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.row >= self.category.rows {
            return None;
        }
        let row = self.row;
        self.row += 1;
        Some(ModelRow {
            category: self.category,
            row,
        })
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.category.rows.saturating_sub(self.row);
        (remaining, Some(remaining))
    }
}

impl ExactSizeIterator for ModelRows<'_> {}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
enum ScalarKey {
    Inapplicable,
    Unknown,
    Text(Arc<str>),
    Integer(i64),
    Float(u64),
}

impl ScalarKey {
    fn new(value: &CifValue) -> Self {
        match value {
            CifValue::Inapplicable => Self::Inapplicable,
            CifValue::Unknown => Self::Unknown,
            CifValue::Text(value) => Self::Text(Arc::clone(value)),
            CifValue::Integer(value) => Self::Integer(*value),
            CifValue::Float(value) => Self::Float(value.to_bits()),
        }
    }
}

fn value_or_unknown<'a>(category: &'a Category, item: &str, row: usize) -> &'a CifValue {
    static UNKNOWN: CifValue = CifValue::Unknown;
    match category.value(item, row) {
        Some(value) => value,
        None => &UNKNOWN,
    }
}
