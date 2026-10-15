//! Streaming category builders that collect CIF cells into compact columns under a memory ceiling.
//!
//! Complexity: amortised O(1) time per cell with hashed category and column lookup;
//! O(cells) space, checked against the ceiling on every push.

use super::column::{ColumnBuilder, PushMemory};
use super::packed::table_bytes;
use super::pairwise::PairwiseStream;
use crate::dense::PAIRWISE;
use crate::{CompactColumn, ModelCategory, ModelCif, ModelCifError, identifier_item};
use molframe_cif::CifScalar;
use std::collections::HashMap;
use std::mem::size_of;

/// Where the streaming dense PAE reader stands.
enum Stream {
    Idle,
    Active(Box<PairwiseStream>),
    Off,
}

pub(super) struct ModelBuilder {
    stream: Stream,
    categories: Vec<CategoryBuilder>,
    lookup: HashMap<Box<str>, usize>,
    last_category: Option<usize>,
    live_bytes: usize,
    limit: usize,
}

impl ModelBuilder {
    pub(super) fn new(limit: usize) -> Self {
        Self {
            stream: Stream::Idle,
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
        if category == PAIRWISE && !matches!(self.stream, Stream::Off) {
            if matches!(self.stream, Stream::Idle) {
                self.stream = Stream::Active(Box::new(PairwiseStream::new(self.categories.len())));
            }
            if let Stream::Active(stream) = &mut self.stream
                && stream.feed(item, value)
            {
                return Ok(());
            }
            return self.divert();
        }
        self.push_compact(category, item, value)
    }

    /// Hands the rows a rejected stream consumed back to the compact table.
    fn divert(&mut self) -> Result<(), ModelCifError> {
        let Stream::Active(stream) = std::mem::replace(&mut self.stream, Stream::Off) else {
            return Ok(());
        };
        self.insert_category(stream.slot(), PAIRWISE)?;
        stream.replay(&mut |item, value| self.push_compact(PAIRWISE, item, value))
    }

    fn insert_category(&mut self, slot: usize, name: &str) -> Result<(), ModelCifError> {
        self.categories
            .try_reserve(1)
            .map_err(|_| ModelCifError::Allocation)?;
        self.lookup
            .try_reserve(1)
            .map_err(|_| ModelCifError::Allocation)?;
        self.categories.insert(slot, CategoryBuilder::new(name));
        for index in self.lookup.values_mut().filter(|index| **index >= slot) {
            *index += 1;
        }
        self.lookup.insert(name.into(), slot);
        self.last_category = None;
        self.refresh_live();
        Ok(())
    }

    fn push_compact(
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
        let mut stream = None;
        if let Stream::Active(mut active) = std::mem::replace(&mut self.stream, Stream::Off) {
            if active.complete() {
                stream = Some(active);
            } else {
                self.stream = Stream::Active(active);
                self.divert()?;
            }
        }
        let limit = self.limit;
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
        let mut model = ModelCif {
            categories,
            pae: None,
        };
        if let Some(stream) = stream {
            attach_stream(&mut model, *stream, limit)?;
        }
        Ok(model)
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
        self.rows = maximum_rows(&self.columns);
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

fn maximum_rows(columns: &[ColumnBuilder]) -> usize {
    let Some(rows) = columns.iter().map(ColumnBuilder::len).max() else {
        return 0;
    };
    rows
}

/// Installs the streamed matrix, or restores the compact table for a non-PAE metric.
fn attach_stream(
    model: &mut ModelCif,
    stream: PairwiseStream,
    limit: usize,
) -> Result<(), ModelCifError> {
    let is_pae = stream.metric_id().is_some_and(|metric| {
        model.confidence().definitions().any(|definition| {
            definition.id == metric.as_str() && definition.metric_type.eq_ignore_ascii_case("PAE")
        })
    });
    if is_pae {
        model.pae = stream.into_dense();
        return Ok(());
    }
    let mut table = CategoryBuilder::new(PAIRWISE);
    stream.replay(&mut |item, value| table.push(item, value).map(|_| ()))?;
    table.normalize()?;
    if table.live_bytes() > limit {
        return Err(ModelCifError::MemoryLimit {
            required: table.live_bytes(),
            limit,
        });
    }
    let slot = stream.slot().min(model.categories.len());
    model.categories.insert(slot, table.finish()?);
    Ok(())
}
