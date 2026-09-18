use super::column::{ColumnBuilder, PushMemory};
use crate::{CompactColumn, ModelCategory, ModelCif, ModelCifError, identifier_item};
use molframe_cif::CifScalar;
use std::collections::HashMap;
use std::mem::size_of;

pub(super) struct ModelBuilder {
    categories: Vec<CategoryBuilder>,
    lookup: HashMap<Box<str>, usize>,
    last_category: Option<usize>,
    live_bytes: usize,
    limit: usize,
}

impl ModelBuilder {
    pub(super) fn new(limit: usize) -> Self {
        Self {
            categories: Vec::new(),
            lookup: HashMap::new(),
            last_category: None,
            live_bytes: 0,
            limit,
        }
    }

    pub(super) fn push(
        &mut self,
        category: &str,
        item: &str,
        value: CifScalar<'_>,
    ) -> Result<(), ModelCifError> {
        let mut peak = self.live_bytes;
        let index = if let Some(index) = self.cached_category(category) {
            index
        } else {
            let before = self.live_bytes;
            let index = self.add_category(category)?;
            self.refresh_live();
            peak = peak.max(before.saturating_add(self.live_bytes));
            index
        };
        let Some(target) = self.categories.get_mut(index) else {
            return Err(ModelCifError::Capacity);
        };
        self.last_category = Some(index);
        let before = target.live_bytes();
        let base = self.live_bytes.saturating_sub(before);
        let memory = target.push(item, value)?;
        self.live_bytes = base.saturating_add(memory.retained);
        peak = peak.max(base.saturating_add(memory.peak));
        self.enforce(peak)
    }

    pub(super) fn finish(mut self) -> Result<ModelCif, ModelCifError> {
        for index in 0..self.categories.len() {
            let before = self.categories[index].live_bytes();
            let base = self.live_bytes.saturating_sub(before);
            let memory = self.categories[index].normalize()?;
            self.live_bytes = base.saturating_add(memory.retained);
            self.enforce(base.saturating_add(memory.peak))?;
        }
        let output_bytes = self
            .categories
            .len()
            .saturating_mul(size_of::<ModelCategory>())
            .saturating_add(
                self.categories
                    .iter()
                    .map(|category| {
                        category
                            .columns
                            .len()
                            .saturating_mul(size_of::<CompactColumn>())
                    })
                    .sum::<usize>(),
            );
        self.enforce(self.live_bytes.saturating_add(output_bytes))?;
        let mut categories = Vec::new();
        categories
            .try_reserve_exact(self.categories.len())
            .map_err(|_| ModelCifError::Allocation)?;
        for category in self.categories {
            categories.push(category.finish()?);
        }
        Ok(ModelCif { categories })
    }

    #[cfg(test)]
    pub(super) const fn tracked_bytes(&self) -> usize {
        self.live_bytes
    }

    fn add_category(&mut self, name: &str) -> Result<usize, ModelCifError> {
        self.categories
            .try_reserve(1)
            .map_err(|_| ModelCifError::Allocation)?;
        self.lookup
            .try_reserve(1)
            .map_err(|_| ModelCifError::Allocation)?;
        let index = self.categories.len();
        self.categories.push(CategoryBuilder::new(name));
        self.lookup.insert(name.into(), index);
        Ok(index)
    }

    fn cached_category(&self, name: &str) -> Option<usize> {
        self.last_category
            .and_then(|index| {
                self.categories
                    .get(index)
                    .filter(|category| category.name.as_ref() == name)
                    .map(|_| index)
            })
            .or_else(|| self.lookup.get(name).copied())
    }

    fn refresh_live(&mut self) {
        let categories = self
            .categories
            .capacity()
            .saturating_mul(size_of::<CategoryBuilder>());
        let names = self.lookup.keys().map(|name| name.len()).sum::<usize>();
        self.live_bytes = categories
            .saturating_add(table_bytes(&self.lookup))
            .saturating_add(names)
            .saturating_add(
                self.categories
                    .iter()
                    .map(CategoryBuilder::live_bytes)
                    .sum::<usize>(),
            );
    }

    fn enforce(&self, required: usize) -> Result<(), ModelCifError> {
        if required > self.limit {
            return Err(ModelCifError::MemoryLimit {
                required,
                limit: self.limit,
            });
        }
        Ok(())
    }
}

struct CategoryBuilder {
    name: Box<str>,
    items: Vec<Box<str>>,
    columns: Vec<ColumnBuilder>,
    lookup: HashMap<Box<str>, usize>,
    last_column: Option<usize>,
    expected_column: usize,
    live_bytes: usize,
    rows: usize,
}

impl CategoryBuilder {
    fn new(name: &str) -> Self {
        Self {
            name: name.into(),
            items: Vec::new(),
            columns: Vec::new(),
            lookup: HashMap::new(),
            last_column: None,
            expected_column: 0,
            live_bytes: name.len(),
            rows: 0,
        }
    }

    const fn live_bytes(&self) -> usize {
        self.live_bytes
    }

    fn push(&mut self, item: &str, value: CifScalar<'_>) -> Result<PushMemory, ModelCifError> {
        let mut peak = self.live_bytes;
        let index = if let Some(index) = self.cached_column(item) {
            index
        } else {
            let before = self.live_bytes;
            let index = self.add_column(item)?;
            self.refresh_live();
            peak = peak.max(before.saturating_add(self.live_bytes));
            index
        };
        let column_count = self.columns.len();
        let Some(column) = self.columns.get_mut(index) else {
            return Err(ModelCifError::Capacity);
        };
        self.last_column = Some(index);
        self.expected_column = (index + 1) % column_count;
        let before = column.live_bytes();
        let base = self.live_bytes.saturating_sub(before);
        let memory = column.push(value)?;
        self.live_bytes = base.saturating_add(memory.retained);
        peak = peak.max(base.saturating_add(memory.peak));
        Ok(PushMemory {
            retained: self.live_bytes,
            peak,
        })
    }

    fn normalize(&mut self) -> Result<PushMemory, ModelCifError> {
        self.rows = match self.columns.iter().map(ColumnBuilder::len).max() {
            Some(rows) => rows,
            None => 0,
        };
        let mut peak = self.live_bytes;
        for index in 0..self.columns.len() {
            while self.columns[index].len() < self.rows {
                let before = self.columns[index].live_bytes();
                let base = self.live_bytes.saturating_sub(before);
                let memory = self.columns[index].push(CifScalar::Unknown)?;
                self.live_bytes = base.saturating_add(memory.retained);
                peak = peak.max(base.saturating_add(memory.peak));
            }
        }
        Ok(PushMemory {
            retained: self.live_bytes,
            peak,
        })
    }

    fn finish(self) -> Result<ModelCategory, ModelCifError> {
        let mut columns = Vec::new();
        columns
            .try_reserve_exact(self.columns.len())
            .map_err(|_| ModelCifError::Allocation)?;
        for column in self.columns {
            let column: CompactColumn = column.finish();
            columns.push(column);
        }
        Ok(ModelCategory {
            name: self.name,
            items: self.items,
            columns,
            rows: self.rows,
        })
    }

    fn add_column(&mut self, item: &str) -> Result<usize, ModelCifError> {
        self.items
            .try_reserve(1)
            .map_err(|_| ModelCifError::Allocation)?;
        self.columns
            .try_reserve(1)
            .map_err(|_| ModelCifError::Allocation)?;
        self.lookup
            .try_reserve(1)
            .map_err(|_| ModelCifError::Allocation)?;
        let index = self.columns.len();
        self.items.push(item.into());
        self.columns
            .push(ColumnBuilder::new(identifier_item(&self.name, item)));
        self.lookup.insert(item.into(), index);
        Ok(index)
    }

    fn cached_column(&self, item: &str) -> Option<usize> {
        self.last_column
            .and_then(|index| self.matching_column(index, item))
            .or_else(|| self.matching_column(self.expected_column, item))
            .or_else(|| self.lookup.get(item).copied())
    }

    fn matching_column(&self, index: usize, item: &str) -> Option<usize> {
        self.items
            .get(index)
            .filter(|candidate| candidate.as_ref() == item)
            .map(|_| index)
    }

    fn refresh_live(&mut self) {
        let item_buffers = self.items.capacity().saturating_mul(size_of::<Box<str>>());
        let column_buffers = self
            .columns
            .capacity()
            .saturating_mul(size_of::<ColumnBuilder>());
        let names = self
            .name
            .len()
            .saturating_add(self.items.iter().map(|item| item.len()).sum::<usize>())
            .saturating_add(self.lookup.keys().map(|item| item.len()).sum::<usize>());
        self.live_bytes = item_buffers
            .saturating_add(column_buffers)
            .saturating_add(table_bytes(&self.lookup))
            .saturating_add(names)
            .saturating_add(
                self.columns
                    .iter()
                    .map(ColumnBuilder::live_bytes)
                    .sum::<usize>(),
            );
    }
}

fn table_bytes<K, V>(values: &HashMap<K, V>) -> usize {
    values.capacity().saturating_mul(
        size_of::<(K, V)>()
            .saturating_add(size_of::<usize>())
            .saturating_add(1),
    )
}
