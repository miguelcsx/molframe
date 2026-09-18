//! Reusable decoded-column storage retained by one `BinaryCIF` source.

use molframe_core::{Code, Diagnostic};

#[derive(Debug)]
pub(super) enum ColumnChunk {
    Integer(Vec<i64>),
    Float(Vec<f32>),
    Text(Vec<u32>),
}

#[derive(Debug)]
pub(super) struct DecodedChunk {
    pub(super) values: ColumnChunk,
    pub(super) mask: Option<Vec<u8>>,
}

impl DecodedChunk {
    pub(super) fn capacity_bytes(&self) -> usize {
        let values = match &self.values {
            ColumnChunk::Integer(values) => capacity_bytes::<i64>(values.capacity()),
            ColumnChunk::Float(values) => capacity_bytes::<f32>(values.capacity()),
            ColumnChunk::Text(values) => capacity_bytes::<u32>(values.capacity()),
        };
        values.saturating_add(
            self.mask
                .as_ref()
                .map_or(0, |mask| capacity_bytes::<u8>(mask.capacity())),
        )
    }

    pub(super) fn additional_capacity_bytes(&self, rows: usize) -> usize {
        let values = match &self.values {
            ColumnChunk::Integer(values) => missing_bytes::<i64>(values.capacity(), rows),
            ColumnChunk::Float(values) => missing_bytes::<f32>(values.capacity(), rows),
            ColumnChunk::Text(values) => missing_bytes::<u32>(values.capacity(), rows),
        };
        values.saturating_add(
            self.mask
                .as_ref()
                .map_or(0, |mask| missing_bytes::<u8>(mask.capacity(), rows)),
        )
    }

    pub(super) fn ensure_capacity(&mut self, rows: usize) -> Result<(), Diagnostic> {
        match &mut self.values {
            ColumnChunk::Integer(values) => reserve(values, rows)?,
            ColumnChunk::Float(values) => reserve(values, rows)?,
            ColumnChunk::Text(values) => reserve(values, rows)?,
        }
        if let Some(mask) = &mut self.mask {
            reserve(mask, rows)?;
        }
        Ok(())
    }

    pub(super) fn clear(&mut self) {
        match &mut self.values {
            ColumnChunk::Integer(values) => values.clear(),
            ColumnChunk::Float(values) => values.clear(),
            ColumnChunk::Text(values) => values.clear(),
        }
        if let Some(mask) = &mut self.mask {
            mask.clear();
        }
    }
}

fn missing_bytes<T>(capacity: usize, rows: usize) -> usize {
    rows.saturating_sub(capacity)
        .saturating_mul(std::mem::size_of::<T>())
}

fn capacity_bytes<T>(capacity: usize) -> usize {
    capacity.saturating_mul(std::mem::size_of::<T>())
}

fn reserve<T>(values: &mut Vec<T>, rows: usize) -> Result<(), Diagnostic> {
    if rows <= values.capacity() {
        return Ok(());
    }
    values
        .try_reserve_exact(rows - values.capacity())
        .map_err(|error| {
            Diagnostic::new(Code::E1902)
                .with_message("BinaryCIF decode workspace allocation failed")
                .with_context("reason", error.to_string())
        })
}
