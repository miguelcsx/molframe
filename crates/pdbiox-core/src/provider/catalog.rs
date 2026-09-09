//! Compact dataset metadata without per-chunk materialisation.

use super::{ChunkDescriptor, ChunkId, DatasetId, LogicalRow, ProviderError};
use std::sync::Arc;

/// Payload carried by a dataset.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum PayloadKind {
    /// Atom topology and base coordinates.
    Structure,
    /// One typed property column.
    Property,
    /// Coordinates for one frame.
    Frame,
    /// Chemical bond rows with globally addressable atom endpoints.
    BondTopology,
}

/// How chunk boundaries are determined.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ChunkLayout {
    /// Every chunk except the tail has the same row count.
    Regular {
        /// Rows in a full chunk.
        rows_per_chunk: u32,
    },
    /// The source provides exact boundaries, while this value guides demand.
    SourceDefined {
        /// Preferred maximum rows requested from the source.
        target_rows: u32,
    },
}

impl ChunkLayout {
    /// Preferred row count used for bounded pull demand.
    #[must_use]
    pub const fn target_rows(self) -> u32 {
        match self {
            Self::Regular { rows_per_chunk } => rows_per_chunk,
            Self::SourceDefined { target_rows } => target_rows,
        }
    }
}

/// O(1) metadata for one logical dataset.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DatasetDescriptor {
    id: DatasetId,
    payload: PayloadKind,
    logical_rows: u64,
    chunk_count: u64,
    first_chunk: ChunkId,
    layout: ChunkLayout,
}

impl DatasetDescriptor {
    /// Describes a regularly chunked logical dataset.
    ///
    /// # Errors
    ///
    /// Returns a typed error for zero chunk size or 64-bit arithmetic overflow.
    pub fn regular(
        id: DatasetId,
        payload: PayloadKind,
        logical_rows: u64,
        rows_per_chunk: u32,
        first_chunk: ChunkId,
    ) -> Result<Self, ProviderError> {
        if rows_per_chunk == 0 {
            return Err(ProviderError::ZeroChunkRows { dataset: id });
        }
        let divisor = u64::from(rows_per_chunk);
        let chunk_count = match logical_rows.checked_add(divisor - 1) {
            Some(rounded) => rounded / divisor,
            None => {
                return Err(ProviderError::IdentityOverflow {
                    dataset: id,
                    field: "chunk count",
                });
            }
        };
        Self::new(
            id,
            payload,
            logical_rows,
            chunk_count,
            first_chunk,
            ChunkLayout::Regular { rows_per_chunk },
        )
    }

    /// Describes source-defined chunk boundaries without storing each boundary.
    ///
    /// # Errors
    ///
    /// Returns a typed error for contradictory counts, zero target size, or
    /// chunk-identity overflow.
    pub fn source_defined(
        id: DatasetId,
        payload: PayloadKind,
        logical_rows: u64,
        chunk_count: u64,
        first_chunk: ChunkId,
        target_rows: u32,
    ) -> Result<Self, ProviderError> {
        Self::new(
            id,
            payload,
            logical_rows,
            chunk_count,
            first_chunk,
            ChunkLayout::SourceDefined { target_rows },
        )
    }

    fn new(
        id: DatasetId,
        payload: PayloadKind,
        logical_rows: u64,
        chunk_count: u64,
        first_chunk: ChunkId,
        layout: ChunkLayout,
    ) -> Result<Self, ProviderError> {
        if layout.target_rows() == 0 {
            return Err(ProviderError::ZeroChunkRows { dataset: id });
        }
        if logical_rows > 0 && chunk_count == 0 {
            return Err(ProviderError::MissingChunks {
                dataset: id,
                rows: logical_rows,
            });
        }
        validate_chunk_range(id, first_chunk, chunk_count)?;
        Ok(Self {
            id,
            payload,
            logical_rows,
            chunk_count,
            first_chunk,
            layout,
        })
    }

    /// Stable dataset identity.
    #[must_use]
    pub const fn id(self) -> DatasetId {
        self.id
    }

    /// Typed payload produced for this dataset.
    #[must_use]
    pub const fn payload(self) -> PayloadKind {
        self.payload
    }

    /// Total logical rows, which may exceed `u32::MAX`.
    #[must_use]
    pub const fn logical_rows(self) -> u64 {
        self.logical_rows
    }

    /// Number of chunks, represented without a per-chunk vector.
    #[must_use]
    pub const fn chunk_count(self) -> u64 {
        self.chunk_count
    }

    /// First stable chunk identity.
    #[must_use]
    pub const fn first_chunk(self) -> ChunkId {
        self.first_chunk
    }

    /// Chunk-boundary policy.
    #[must_use]
    pub const fn layout(self) -> ChunkLayout {
        self.layout
    }

    /// Resolves a regular chunk in constant time.
    ///
    /// # Errors
    ///
    /// Returns [`ProviderError::UnknownChunk`] for a foreign chunk or when the
    /// source owns its own boundaries.
    pub fn regular_chunk(self, chunk: ChunkId) -> Result<ChunkDescriptor, ProviderError> {
        let ChunkLayout::Regular { rows_per_chunk } = self.layout else {
            return Err(ProviderError::UnknownChunk {
                dataset: self.id,
                chunk,
            });
        };
        let ordinal = self.chunk_ordinal(chunk)?;
        let start = ordinal.checked_mul(u64::from(rows_per_chunk)).ok_or(
            ProviderError::IdentityOverflow {
                dataset: self.id,
                field: "logical row",
            },
        )?;
        let remaining =
            self.logical_rows
                .checked_sub(start)
                .ok_or(ProviderError::LogicalRowOutOfRange {
                    dataset: self.id,
                    row: LogicalRow::new(start),
                })?;
        let rows = remaining.min(u64::from(rows_per_chunk));
        let rows = u32::try_from(rows).map_err(|_| ProviderError::IdentityOverflow {
            dataset: self.id,
            field: "local row count",
        })?;
        ChunkDescriptor::new(self.id, chunk, LogicalRow::new(start), rows)
    }

    pub(crate) fn chunk_ordinal(self, chunk: ChunkId) -> Result<u64, ProviderError> {
        let Some(ordinal) = chunk.get().checked_sub(self.first_chunk.get()) else {
            return Err(ProviderError::UnknownChunk {
                dataset: self.id,
                chunk,
            });
        };
        if ordinal >= self.chunk_count {
            return Err(ProviderError::UnknownChunk {
                dataset: self.id,
                chunk,
            });
        }
        Ok(ordinal)
    }
}

/// Hierarchical catalog whose memory scales with datasets, not logical chunks.
#[derive(Clone, Debug, Default)]
pub struct DatasetCatalog {
    datasets: Arc<[DatasetDescriptor]>,
}

impl DatasetCatalog {
    /// Validates, sorts, and shares dataset metadata.
    ///
    /// # Errors
    ///
    /// Rejects duplicate dataset identities and overlapping chunk namespaces.
    pub fn new(mut datasets: Vec<DatasetDescriptor>) -> Result<Self, ProviderError> {
        datasets.sort_unstable_by_key(|dataset| dataset.id());
        validate_dataset_ids(&datasets)?;
        validate_chunk_ranges(&datasets)?;
        Ok(Self {
            datasets: datasets.into(),
        })
    }

    /// Number of logical datasets.
    #[must_use]
    pub fn len(&self) -> usize {
        self.datasets.len()
    }

    /// Reports whether the catalog contains no datasets.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.datasets.is_empty()
    }

    /// Shared sorted descriptors.
    #[must_use]
    pub fn datasets(&self) -> &[DatasetDescriptor] {
        &self.datasets
    }

    /// Finds one dataset in logarithmic time.
    #[must_use]
    pub fn get(&self, id: DatasetId) -> Option<DatasetDescriptor> {
        self.datasets
            .binary_search_by_key(&id, |dataset| dataset.id())
            .ok()
            .and_then(|index| self.datasets.get(index).copied())
    }
}

fn validate_chunk_range(
    dataset: DatasetId,
    first: ChunkId,
    count: u64,
) -> Result<(), ProviderError> {
    if count == 0 {
        return Ok(());
    }
    first
        .get()
        .checked_add(count - 1)
        .map(|_| ())
        .ok_or(ProviderError::IdentityOverflow {
            dataset,
            field: "chunk identity",
        })
}

fn validate_dataset_ids(datasets: &[DatasetDescriptor]) -> Result<(), ProviderError> {
    for pair in datasets.windows(2) {
        if pair[0].id() == pair[1].id() {
            return Err(ProviderError::DuplicateDataset {
                dataset: pair[0].id(),
            });
        }
    }
    Ok(())
}

fn validate_chunk_ranges(datasets: &[DatasetDescriptor]) -> Result<(), ProviderError> {
    let mut ranges: Vec<_> = datasets
        .iter()
        .filter(|dataset| dataset.chunk_count() > 0)
        .map(|dataset| {
            let end = dataset.first_chunk().get() + dataset.chunk_count() - 1;
            (dataset.first_chunk().get(), end)
        })
        .collect();
    ranges.sort_unstable_by_key(|range| range.0);
    for pair in ranges.windows(2) {
        if pair[1].0 <= pair[0].1 {
            return Err(ProviderError::OverlappingChunkRange {
                chunk: ChunkId::new(pair[1].0),
            });
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "catalog_tests.rs"]
mod tests;
