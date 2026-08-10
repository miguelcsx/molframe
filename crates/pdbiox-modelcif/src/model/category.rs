use pdbiox_cif::{Category, CifValue};

/// One `ModelCIF` row with dictionary item names retained.
#[derive(Clone, Debug, PartialEq)]
pub struct ModelRow {
    values: Vec<CifValue>,
}

impl ModelRow {
    pub(crate) fn values(&self) -> &[CifValue] {
        &self.values
    }

    /// One value in item order.
    #[must_use]
    pub fn value(&self, index: usize) -> Option<&CifValue> {
        self.values.get(index)
    }
}

/// One complete `ma_*` category, suitable for lossless semantic rewriting.
#[derive(Clone, Debug, PartialEq)]
pub struct ModelCategory {
    name: Box<str>,
    items: Vec<Box<str>>,
    rows: Vec<ModelRow>,
}

impl ModelCategory {
    pub(crate) fn from_category(category: &Category) -> Self {
        let items: Vec<Box<str>> = category.items().map(Into::into).collect();
        let rows = (0..category.row_count())
            .map(|row| ModelRow {
                values: items
                    .iter()
                    .map(|item| match category.value(item, row) {
                        Some(value) => value.clone(),
                        None => CifValue::Unknown,
                    })
                    .collect(),
            })
            .collect();
        Self {
            name: category.name().into(),
            items,
            rows,
        }
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

    /// Rows in source order.
    #[must_use]
    pub fn rows(&self) -> &[ModelRow] {
        &self.rows
    }

    /// Looks up a value by dictionary item and row.
    #[must_use]
    pub fn value(&self, item: &str, row: usize) -> Option<&CifValue> {
        let column = self
            .items
            .iter()
            .position(|candidate| candidate.as_ref() == item)?;
        self.rows.get(row)?.value(column)
    }
}
