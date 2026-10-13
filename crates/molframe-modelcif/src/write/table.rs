//! Cell-addressed view of the retained `ma_*` tables for binary writers.
//!
//! Compact categories and the dense PAE matrix are exposed as tables whose
//! cells are borrowed on demand, so a writer can encode columns straight from
//! storage without rendering text or rebuilding rows.

use crate::{DensePae, ModelCategory, ModelCif};

/// One borrowed cell of an [`ExportTable`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ExportCell<'a> {
    /// The CIF `.` value.
    Inapplicable,
    /// The CIF `?` value.
    Unknown,
    /// Text.
    Text(&'a str),
    /// Integer.
    Integer(i64),
    /// Floating point.
    Float(f64),
}

#[derive(Clone, Copy, Debug)]
enum Role {
    Ordinal,
    Model,
    Metric,
    Chain(bool),
    Sequence(bool),
    Component(bool),
    Value,
}

#[derive(Debug)]
enum Source<'a> {
    Compact(&'a ModelCategory),
    Dense { pae: &'a DensePae, roles: Vec<Role> },
}

/// A retained table addressed by column and row.
#[derive(Debug)]
pub struct ExportTable<'a> {
    source: Source<'a>,
}

impl<'a> ExportTable<'a> {
    fn compact(category: &'a ModelCategory) -> Self {
        Self {
            source: Source::Compact(category),
        }
    }

    fn dense(pae: &'a DensePae) -> Self {
        let roles = pae
            .items()
            .iter()
            .map(|item| match &**item {
                "ordinal_id" => Role::Ordinal,
                "model_id" => Role::Model,
                "metric_id" => Role::Metric,
                "label_asym_id_1" => Role::Chain(false),
                "label_asym_id_2" => Role::Chain(true),
                "label_seq_id_1" => Role::Sequence(false),
                "label_seq_id_2" => Role::Sequence(true),
                "label_comp_id_1" => Role::Component(false),
                "label_comp_id_2" => Role::Component(true),
                _ => Role::Value,
            })
            .collect();
        Self {
            source: Source::Dense { pae, roles },
        }
    }

    /// Category name without a leading underscore.
    #[must_use]
    pub fn name(&self) -> &str {
        match &self.source {
            Source::Compact(category) => category.name(),
            Source::Dense { .. } => "ma_qa_metric_local_pairwise",
        }
    }

    /// Number of rows.
    #[must_use]
    pub fn row_count(&self) -> usize {
        match &self.source {
            Source::Compact(category) => category.row_count(),
            Source::Dense { pae, .. } => {
                let size = pae.matrix().size();
                size.saturating_mul(size)
            }
        }
    }

    /// Number of columns.
    #[must_use]
    pub fn column_count(&self) -> usize {
        match &self.source {
            Source::Compact(category) => category.items().len(),
            Source::Dense { roles, .. } => roles.len(),
        }
    }

    /// Item name of one column, empty when out of range.
    #[must_use]
    pub fn item(&self, column: usize) -> &str {
        let items = match &self.source {
            Source::Compact(category) => category.items(),
            Source::Dense { pae, .. } => pae.items(),
        };
        items.get(column).map_or("", |item| item)
    }

    /// One cell; out-of-range cells read as unknown.
    #[must_use]
    pub fn cell(&self, column: usize, row: usize) -> ExportCell<'_> {
        match &self.source {
            Source::Compact(category) => compact_cell(category, column, row),
            Source::Dense { pae, roles } => roles
                .get(column)
                .map_or(ExportCell::Unknown, |role| dense_cell(pae, *role, row)),
        }
    }
}

fn compact_cell(category: &ModelCategory, column: usize, row: usize) -> ExportCell<'_> {
    use crate::ValueRef;
    match category.value_ref(column, row) {
        Some(ValueRef::Inapplicable) => ExportCell::Inapplicable,
        Some(ValueRef::Text(value)) => ExportCell::Text(value),
        Some(ValueRef::Integer(value)) => ExportCell::Integer(value),
        Some(ValueRef::Float(value)) => ExportCell::Float(value),
        Some(ValueRef::Unknown) | None => ExportCell::Unknown,
    }
}

fn dense_cell(pae: &DensePae, role: Role, pair: usize) -> ExportCell<'_> {
    let size = pae.matrix().size();
    if size == 0 {
        return ExportCell::Unknown;
    }
    let (row, column) = (pair / size, pair % size);
    let pick = |second: bool| if second { column } else { row };
    match role {
        Role::Ordinal => i64::try_from(pair + 1).map_or(ExportCell::Unknown, ExportCell::Integer),
        Role::Model => ExportCell::Text(pae.model_id()),
        Role::Metric => ExportCell::Text(pae.metric_id()),
        Role::Chain(second) => pae
            .residues()
            .get(pick(second))
            .map_or(ExportCell::Unknown, |residue| {
                ExportCell::Text(residue.chain)
            }),
        Role::Sequence(second) => pae
            .residues()
            .get(pick(second))
            .map_or(ExportCell::Unknown, |residue| {
                ExportCell::Integer(residue.sequence)
            }),
        Role::Component(second) => pae
            .residues()
            .get(pick(second))
            .and_then(|residue| residue.component)
            .map_or(ExportCell::Unknown, ExportCell::Text),
        Role::Value => pae
            .matrix()
            .get_f64(row, column)
            .map_or(ExportCell::Unknown, ExportCell::Float),
    }
}

impl ModelCif {
    /// Retained tables in source order, with the dense PAE where its table sat.
    ///
    /// Empty compact categories are skipped. Nothing is copied: each table
    /// borrows the stored columns or matrix.
    #[must_use]
    pub fn export_tables(&self) -> Vec<ExportTable<'_>> {
        let mut tables = Vec::with_capacity(self.categories.len() + 1);
        let dense_at = |index: usize| {
            self.dense_pae()
                .filter(|dense| dense.slot == index && dense.matrix().size() > 0)
                .map(ExportTable::dense)
        };
        for (index, category) in self.categories.iter().enumerate() {
            tables.extend(dense_at(index));
            if category.row_count() > 0 {
                tables.push(ExportTable::compact(category));
            }
        }
        tables.extend(dense_at(self.categories.len()));
        tables
    }
}

#[cfg(test)]
#[path = "table_tests.rs"]
mod tests;
