//! Allocation-stable transactional checkpoints for streaming column decoders.

use super::{ColumnStream, Expansion, FloatInner, FloatStream, IntegerStream, ValueStream};

#[derive(Debug)]
pub(in crate::reader::batch) struct ColumnCheckpoints {
    columns: Vec<ColumnCheckpoint>,
    expansions: Vec<Expansion>,
}

#[derive(Clone, Copy, Debug)]
struct StreamCheckpoint {
    offset: u64,
    expansion_start: usize,
    expansion_end: usize,
}

#[derive(Clone, Copy, Debug)]
struct ColumnCheckpoint {
    values: StreamCheckpoint,
    mask: Option<StreamCheckpoint>,
}

impl ColumnCheckpoints {
    pub(in crate::reader::batch) fn new(columns: &[ColumnStream]) -> Self {
        let expansion_count = columns
            .iter()
            .map(ColumnStream::expansion_count)
            .sum::<usize>();
        Self {
            columns: Vec::with_capacity(columns.len()),
            expansions: Vec::with_capacity(expansion_count),
        }
    }

    pub(in crate::reader::batch) fn capture(&mut self, columns: &[ColumnStream]) {
        self.columns.clear();
        self.expansions.clear();
        for column in columns {
            self.columns.push(column.checkpoint(&mut self.expansions));
        }
    }

    pub(in crate::reader::batch) fn restore(&self, columns: &mut [ColumnStream]) {
        for (column, checkpoint) in columns.iter_mut().zip(&self.columns) {
            column.restore(*checkpoint, &self.expansions);
        }
    }

    pub(in crate::reader::batch) fn retained_bytes(&self) -> usize {
        self.columns
            .capacity()
            .saturating_mul(std::mem::size_of::<ColumnCheckpoint>())
            .saturating_add(
                self.expansions
                    .capacity()
                    .saturating_mul(std::mem::size_of::<Expansion>()),
            )
    }
}

impl ColumnStream {
    fn expansion_count(&self) -> usize {
        self.values
            .expansion_count()
            .saturating_add(self.mask.as_ref().map_or(0, IntegerStream::expansion_count))
    }

    fn checkpoint(&self, expansions: &mut Vec<Expansion>) -> ColumnCheckpoint {
        ColumnCheckpoint {
            values: self.values.checkpoint(expansions),
            mask: self
                .mask
                .as_ref()
                .map(|stream| stream.checkpoint(expansions)),
        }
    }

    fn restore(&mut self, checkpoint: ColumnCheckpoint, expansions: &[Expansion]) {
        self.values.restore(checkpoint.values, expansions);
        if let (Some(mask), Some(saved)) = (&mut self.mask, checkpoint.mask) {
            mask.restore(saved, expansions);
        }
    }
}

impl ValueStream {
    fn expansion_count(&self) -> usize {
        match self {
            Self::Integer(stream) => stream.expansion_count(),
            Self::Float(stream) => stream.expansion_count(),
            Self::Text(stream) => stream.indices.expansion_count(),
        }
    }

    fn checkpoint(&self, expansions: &mut Vec<Expansion>) -> StreamCheckpoint {
        match self {
            Self::Integer(stream) => stream.checkpoint(expansions),
            Self::Float(stream) => stream.checkpoint(expansions),
            Self::Text(stream) => stream.indices.checkpoint(expansions),
        }
    }

    fn restore(&mut self, checkpoint: StreamCheckpoint, expansions: &[Expansion]) {
        match self {
            Self::Integer(stream) => stream.restore(checkpoint, expansions),
            Self::Float(stream) => stream.restore(checkpoint, expansions),
            Self::Text(stream) => stream.indices.restore(checkpoint, expansions),
        }
    }
}

impl IntegerStream {
    fn expansion_count(&self) -> usize {
        self.expansions.len()
    }

    fn checkpoint(&self, expansions: &mut Vec<Expansion>) -> StreamCheckpoint {
        let expansion_start = expansions.len();
        expansions.extend_from_slice(&self.expansions);
        StreamCheckpoint {
            offset: self.input.checkpoint(),
            expansion_start,
            expansion_end: expansions.len(),
        }
    }

    fn restore(&mut self, checkpoint: StreamCheckpoint, expansions: &[Expansion]) {
        self.input.restore(checkpoint.offset);
        self.expansions.clear();
        if let Some(saved) = expansions.get(checkpoint.expansion_start..checkpoint.expansion_end) {
            self.expansions.extend_from_slice(saved);
        }
    }
}

impl FloatStream {
    fn expansion_count(&self) -> usize {
        match &self.inner {
            FloatInner::Literal { .. } => 0,
            FloatInner::Integer(stream) => stream.expansion_count(),
        }
    }

    fn checkpoint(&self, expansions: &mut Vec<Expansion>) -> StreamCheckpoint {
        match &self.inner {
            FloatInner::Literal { input, .. } => StreamCheckpoint {
                offset: input.checkpoint(),
                expansion_start: expansions.len(),
                expansion_end: expansions.len(),
            },
            FloatInner::Integer(stream) => stream.checkpoint(expansions),
        }
    }

    fn restore(&mut self, checkpoint: StreamCheckpoint, expansions: &[Expansion]) {
        match &mut self.inner {
            FloatInner::Literal { input, .. } => input.restore(checkpoint.offset),
            FloatInner::Integer(stream) => stream.restore(checkpoint, expansions),
        }
    }
}
