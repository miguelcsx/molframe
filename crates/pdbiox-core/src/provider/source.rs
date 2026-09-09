//! Declarative adapters from immutable structure snapshots to provider chunks.

use super::{
    BondChunk, ChunkDescriptor, ChunkId, DatasetDescriptor, DatasetId, FrameChunk, LogicalRow,
    PayloadKind, PropertyChunk, ProviderError, StructureChunk,
};
use crate::chunk::TARGET_CHUNK_ATOMS;
use crate::structure::Structure;
use std::sync::Arc;

/// Default maximum number of bonds independently requested in one payload.
pub const TARGET_CHUNK_BONDS: u32 = 65_536;

/// Zero-copy provider for regularly paged chemical bond topology.
#[derive(Clone, Debug)]
pub struct BondChunkProvider {
    structure: Structure,
    dataset: DatasetDescriptor,
    atom_dataset: DatasetId,
    atom_logical_start: LogicalRow,
}

impl BondChunkProvider {
    /// Creates a provider using the default bounded bond payload size.
    ///
    /// # Errors
    ///
    /// Rejects unavailable topology and identity or chunk-range overflow.
    pub fn new(
        dataset: DatasetId,
        first_chunk: ChunkId,
        atom_dataset: DatasetId,
        atom_logical_start: LogicalRow,
        structure: Structure,
    ) -> Result<Self, ProviderError> {
        Self::with_rows_per_chunk(
            dataset,
            first_chunk,
            atom_dataset,
            atom_logical_start,
            structure,
            TARGET_CHUNK_BONDS,
        )
    }

    /// Creates a provider with an explicit bounded payload row count.
    ///
    /// # Errors
    ///
    /// Rejects zero chunk rows, unavailable topology, or identity overflow.
    pub fn with_rows_per_chunk(
        dataset: DatasetId,
        first_chunk: ChunkId,
        atom_dataset: DatasetId,
        atom_logical_start: LogicalRow,
        structure: Structure,
        rows_per_chunk: u32,
    ) -> Result<Self, ProviderError> {
        if !structure.data().bonds.is_available() {
            return Err(ProviderError::BondTopologyUnavailable);
        }
        validate_atom_endpoint_range(atom_dataset, atom_logical_start, structure.atom_count())?;
        let logical_rows = u64::try_from(structure.data().bonds.len()).map_err(|_| {
            ProviderError::IdentityOverflow {
                dataset,
                field: "bond table rows",
            }
        })?;
        let descriptor = DatasetDescriptor::regular(
            dataset,
            PayloadKind::BondTopology,
            logical_rows,
            rows_per_chunk,
            first_chunk,
        )?;
        Ok(Self {
            structure,
            dataset: descriptor,
            atom_dataset,
            atom_logical_start,
        })
    }

    /// Compact dataset metadata; no per-chunk index is retained.
    #[must_use]
    pub const fn dataset(&self) -> DatasetDescriptor {
        self.dataset
    }

    /// Produces one shared topology payload in O(1), without a table scan.
    ///
    /// # Errors
    ///
    /// Rejects chunk identities outside this provider or source index overflow.
    pub fn chunk(&self, chunk: ChunkId) -> Result<BondChunk, ProviderError> {
        let descriptor = self.dataset.regular_chunk(chunk)?;
        let start = u32::try_from(descriptor.logical_start().get()).map_err(|_| {
            ProviderError::IdentityOverflow {
                dataset: self.dataset.id(),
                field: "bond source row",
            }
        })?;
        let end = start
            .checked_add(descriptor.rows())
            .ok_or(ProviderError::IdentityOverflow {
                dataset: self.dataset.id(),
                field: "bond source range",
            })?;
        BondChunk::shared(
            descriptor,
            self.structure.clone(),
            start..end,
            self.atom_dataset,
            self.atom_logical_start,
        )
    }
}

fn validate_atom_endpoint_range(
    dataset: DatasetId,
    first: LogicalRow,
    rows: u32,
) -> Result<(), ProviderError> {
    if rows == 0 {
        return Ok(());
    }
    first
        .get()
        .checked_add(u64::from(rows - 1))
        .map(|_| ())
        .ok_or(ProviderError::IdentityOverflow {
            dataset,
            field: "atom endpoint range",
        })
}

/// Zero-copy provider for structure chunks already present in one snapshot.
#[derive(Clone, Debug)]
pub struct StructureChunkProvider {
    structure: Structure,
    dataset: DatasetDescriptor,
}

impl StructureChunkProvider {
    /// Creates a provider whose chunk namespace begins at `first_chunk`.
    ///
    /// # Errors
    ///
    /// Returns a typed error when chunk count or identity arithmetic overflows.
    pub fn new(
        dataset: DatasetId,
        first_chunk: ChunkId,
        structure: Structure,
    ) -> Result<Self, ProviderError> {
        let chunk_count = u64::try_from(structure.data().chunks.len()).map_err(|_| {
            ProviderError::IdentityOverflow {
                dataset,
                field: "structure chunk count",
            }
        })?;
        let descriptor = DatasetDescriptor::source_defined(
            dataset,
            PayloadKind::Structure,
            u64::from(structure.atom_count()),
            chunk_count,
            first_chunk,
            TARGET_CHUNK_ATOMS,
        )?;
        Ok(Self {
            structure,
            dataset: descriptor,
        })
    }

    /// Compact dataset metadata.
    #[must_use]
    pub const fn dataset(&self) -> DatasetDescriptor {
        self.dataset
    }

    /// Produces one shared structure payload in constant time.
    ///
    /// # Errors
    ///
    /// Rejects chunk identities outside this provider.
    pub fn chunk(&self, chunk: ChunkId) -> Result<StructureChunk, ProviderError> {
        let ordinal = self.dataset.chunk_ordinal(chunk)?;
        let index = usize::try_from(ordinal).map_err(|_| ProviderError::IdentityOverflow {
            dataset: self.dataset.id(),
            field: "storage chunk index",
        })?;
        let storage = self
            .structure
            .data()
            .chunks
            .get(index)
            .ok_or(ProviderError::StorageChunkUnavailable { index })?;
        let descriptor = ChunkDescriptor::new(
            self.dataset.id(),
            chunk,
            super::LogicalRow::new(u64::from(storage.atoms().start)),
            storage.len(),
        )?;
        StructureChunk::shared(descriptor, self.structure.clone(), index)
    }
}

/// Zero-copy provider for one typed atom property.
#[derive(Clone, Debug)]
pub struct PropertyChunkProvider {
    structure: Structure,
    name: Arc<str>,
    dataset: DatasetDescriptor,
}

impl PropertyChunkProvider {
    /// Creates a provider for the canonical atom-aligned pLDDT confidence column.
    ///
    /// # Errors
    ///
    /// Rejects structures without lowered pLDDT, row misalignment, or identity overflow.
    pub fn plddt(
        dataset: DatasetId,
        first_chunk: ChunkId,
        structure: Structure,
    ) -> Result<Self, ProviderError> {
        Self::new(
            dataset,
            first_chunk,
            structure,
            crate::annotation::PLDDT_ANNOTATION,
        )
    }

    /// Creates a property provider over the structure's existing chunk layout.
    ///
    /// # Errors
    ///
    /// Rejects absent properties, row misalignment, or identity overflow.
    pub fn new(
        dataset: DatasetId,
        first_chunk: ChunkId,
        structure: Structure,
        name: impl Into<Arc<str>>,
    ) -> Result<Self, ProviderError> {
        let name = name.into();
        let property = structure.annotations().get(&name).ok_or_else(|| {
            ProviderError::PropertyUnavailable {
                name: Box::from(name.as_ref()),
            }
        })?;
        if property.len() != structure.atom_count() {
            return Err(ProviderError::PropertyLengthMismatch {
                name: Box::from(name.as_ref()),
                expected: structure.atom_count(),
                actual: property.len(),
            });
        }
        let chunk_count = u64::try_from(structure.data().chunks.len()).map_err(|_| {
            ProviderError::IdentityOverflow {
                dataset,
                field: "property chunk count",
            }
        })?;
        let descriptor = DatasetDescriptor::source_defined(
            dataset,
            PayloadKind::Property,
            u64::from(structure.atom_count()),
            chunk_count,
            first_chunk,
            TARGET_CHUNK_ATOMS,
        )?;
        Ok(Self {
            structure,
            name,
            dataset: descriptor,
        })
    }

    /// Compact dataset metadata.
    #[must_use]
    pub const fn dataset(&self) -> DatasetDescriptor {
        self.dataset
    }

    /// Produces one property payload without copying its values.
    ///
    /// # Errors
    ///
    /// Rejects chunk identities outside this provider.
    pub fn chunk(&self, chunk: ChunkId) -> Result<PropertyChunk, ProviderError> {
        let ordinal = self.dataset.chunk_ordinal(chunk)?;
        let index = usize::try_from(ordinal).map_err(|_| ProviderError::IdentityOverflow {
            dataset: self.dataset.id(),
            field: "storage chunk index",
        })?;
        let storage = self
            .structure
            .data()
            .chunks
            .get(index)
            .ok_or(ProviderError::StorageChunkUnavailable { index })?;
        let descriptor = ChunkDescriptor::new(
            self.dataset.id(),
            chunk,
            super::LogicalRow::new(u64::from(storage.atoms().start)),
            storage.len(),
        )?;
        PropertyChunk::shared(
            descriptor,
            self.structure.clone(),
            Arc::clone(&self.name),
            index,
        )
    }
}

/// Zero-copy provider for one coordinate frame over an existing chunk layout.
#[derive(Clone, Debug)]
pub struct FrameChunkProvider {
    structure: Structure,
    model: crate::ModelIndex,
    storage_chunks: std::ops::Range<usize>,
    dataset: DatasetDescriptor,
}

impl FrameChunkProvider {
    /// Creates a frame provider for one directly stored model.
    ///
    /// # Errors
    ///
    /// Rejects ragged or absent model storage and identity overflow.
    pub fn new(
        dataset: DatasetId,
        first_chunk: ChunkId,
        structure: Structure,
        model: crate::ModelIndex,
    ) -> Result<Self, ProviderError> {
        let Some(coordinates) = structure.data().coords.block(model) else {
            return Err(ProviderError::ModelUnavailable { model: model.get() });
        };
        let storage_chunks = model_chunk_range(&structure, model);
        let chunk_count =
            u64::try_from(storage_chunks.len()).map_err(|_| ProviderError::IdentityOverflow {
                dataset,
                field: "frame chunk count",
            })?;
        let descriptor = DatasetDescriptor::source_defined(
            dataset,
            PayloadKind::Frame,
            u64::from(coordinates.len()),
            chunk_count,
            first_chunk,
            TARGET_CHUNK_ATOMS,
        )?;
        Ok(Self {
            structure,
            model,
            storage_chunks,
            dataset: descriptor,
        })
    }

    /// Compact dataset metadata.
    #[must_use]
    pub const fn dataset(&self) -> DatasetDescriptor {
        self.dataset
    }

    /// Produces one coordinate payload sharing the frame's aligned storage.
    ///
    /// # Errors
    ///
    /// Rejects chunk identities outside this provider or inconsistent models.
    pub fn chunk(&self, chunk: ChunkId) -> Result<FrameChunk, ProviderError> {
        let ordinal = self.dataset.chunk_ordinal(chunk)?;
        let local = usize::try_from(ordinal).map_err(|_| ProviderError::IdentityOverflow {
            dataset: self.dataset.id(),
            field: "storage chunk index",
        })?;
        let index = self.storage_chunks.start.checked_add(local).ok_or(
            ProviderError::IdentityOverflow {
                dataset: self.dataset.id(),
                field: "storage chunk index",
            },
        )?;
        if index >= self.storage_chunks.end {
            return Err(ProviderError::UnknownChunk {
                dataset: self.dataset.id(),
                chunk,
            });
        }
        let storage = &self.structure.data().chunks[index];
        let range = storage.atoms();
        let descriptor = ChunkDescriptor::new(
            self.dataset.id(),
            chunk,
            super::LogicalRow::new(u64::from(range.start)),
            storage.len(),
        )?;
        let coordinates = self
            .structure
            .data()
            .coords
            .block(self.model)
            .cloned()
            .ok_or(ProviderError::ModelUnavailable {
                model: self.model.get(),
            })?;
        FrameChunk::shared(descriptor, coordinates, range)
    }
}

fn model_chunk_range(structure: &Structure, model: crate::ModelIndex) -> std::ops::Range<usize> {
    let chunks = structure.data().chunks.as_slice();
    let start = chunks.partition_point(|chunk| chunk.model() < model.get());
    let end = chunks.partition_point(|chunk| chunk.model() <= model.get());
    start..end
}

#[cfg(test)]
#[path = "source_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "bond_source_tests.rs"]
mod bond_tests;
