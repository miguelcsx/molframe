//! Typed chunk payloads with shared ownership and borrowed projections.

use super::{ChunkId, DatasetId, LocalRow, LogicalRow, ProviderError};
use crate::annotation::AtomAnnotation;
use crate::chunk::AtomChunk;
use crate::column::Presence;
use crate::coords::CoordinateBlock;
use crate::structure::Structure;
use std::ops::Range;
use std::sync::Arc;

/// Stable global identity and compact local extent of one chunk.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ChunkDescriptor {
    dataset: DatasetId,
    chunk: ChunkId,
    logical_start: LogicalRow,
    rows: u32,
}

impl ChunkDescriptor {
    /// Creates a chunk descriptor.
    ///
    /// # Errors
    ///
    /// Rejects a logical range whose end exceeds `u64`.
    pub fn new(
        dataset: DatasetId,
        chunk: ChunkId,
        logical_start: LogicalRow,
        rows: u32,
    ) -> Result<Self, ProviderError> {
        if rows > 0
            && logical_start
                .get()
                .checked_add(u64::from(rows - 1))
                .is_none()
        {
            return Err(ProviderError::IdentityOverflow {
                dataset,
                field: "logical row range",
            });
        }
        Ok(Self {
            dataset,
            chunk,
            logical_start,
            rows,
        })
    }

    /// Dataset identity.
    #[must_use]
    pub const fn dataset(self) -> DatasetId {
        self.dataset
    }

    /// Chunk identity.
    #[must_use]
    pub const fn chunk(self) -> ChunkId {
        self.chunk
    }

    /// First global row represented by local row zero.
    #[must_use]
    pub const fn logical_start(self) -> LogicalRow {
        self.logical_start
    }

    /// Number of locally addressable rows.
    #[must_use]
    pub const fn rows(self) -> u32 {
        self.rows
    }

    /// Resolves a compact local row to its stable global row.
    ///
    /// # Errors
    ///
    /// Returns a typed error when the local row lies outside this chunk.
    pub fn resolve(self, local: LocalRow) -> Result<LogicalRow, ProviderError> {
        if local.get() >= self.rows {
            return Err(ProviderError::LocalRowOutOfRange {
                row: local.get(),
                rows: self.rows,
            });
        }
        let value = self
            .logical_start
            .get()
            .checked_add(u64::from(local.get()))
            .ok_or(ProviderError::IdentityOverflow {
                dataset: self.dataset,
                field: "logical row",
            })?;
        Ok(LogicalRow::new(value))
    }
}

/// One structure chunk retaining its immutable native snapshot.
#[derive(Clone, Debug)]
pub struct StructureChunk {
    descriptor: ChunkDescriptor,
    structure: Structure,
    storage_chunk: usize,
    coordinates: CoordinateBlock,
    coordinate_range: Range<u32>,
}

impl StructureChunk {
    /// Adapts one existing storage chunk without copying coordinates.
    ///
    /// # Errors
    ///
    /// Rejects missing chunks, models, coordinate ranges, or row mismatches.
    pub fn shared(
        descriptor: ChunkDescriptor,
        structure: Structure,
        storage_chunk: usize,
    ) -> Result<Self, ProviderError> {
        let chunk = structure.data().chunks.get(storage_chunk).ok_or(
            ProviderError::StorageChunkUnavailable {
                index: storage_chunk,
            },
        )?;
        if chunk.len() != descriptor.rows() {
            return Err(ProviderError::PayloadLengthMismatch {
                declared: descriptor.rows(),
                actual: chunk.len(),
            });
        }
        let model = crate::ModelIndex::new(chunk.model());
        let coordinates = structure.data().coords.block(model).cloned().ok_or(
            ProviderError::ModelUnavailable {
                model: chunk.model(),
            },
        )?;
        let coordinate_range = chunk.atoms();
        if coordinates.range(coordinate_range.clone()).is_none() {
            return Err(ProviderError::StorageChunkUnavailable {
                index: storage_chunk,
            });
        }
        Ok(Self {
            descriptor,
            structure,
            storage_chunk,
            coordinates,
            coordinate_range,
        })
    }

    /// Stable chunk metadata.
    #[must_use]
    pub const fn descriptor(&self) -> ChunkDescriptor {
        self.descriptor
    }

    /// Columnar atom metadata borrowed from the retained snapshot.
    #[must_use]
    pub fn atoms(&self) -> &AtomChunk {
        &self.structure.data().chunks[self.storage_chunk]
    }

    /// Coordinates borrowed directly from shared native storage.
    #[must_use]
    pub fn positions(&self) -> &[[f32; 3]] {
        match self.coordinates.range(self.coordinate_range.clone()) {
            Some(positions) => positions,
            None => &[],
        }
    }

    /// Retained immutable structure snapshot.
    #[must_use]
    pub const fn structure(&self) -> &Structure {
        &self.structure
    }
}

/// Physical type of a property payload.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PropertyKind {
    /// Boolean values.
    Boolean,
    /// Signed 64-bit integers.
    Integer,
    /// IEEE-754 64-bit values.
    Real,
    /// Interned symbols represented by their dictionary identity.
    Symbol,
}

/// One property value projected without materialising a column.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PropertyValue {
    /// Boolean value and validity.
    Boolean(bool, Presence),
    /// Integer value and validity.
    Integer(i64, Presence),
    /// Real value and validity.
    Real(f64, Presence),
    /// Interned symbol and validity.
    Symbol(crate::SymbolId, Presence),
}

/// One typed property chunk retaining its immutable source snapshot.
#[derive(Clone, Debug)]
pub struct PropertyChunk {
    descriptor: ChunkDescriptor,
    structure: Structure,
    name: Arc<str>,
    atom_start: u32,
    storage_chunk: usize,
}

impl PropertyChunk {
    /// Adapts a property over one existing structure chunk.
    ///
    /// # Errors
    ///
    /// Rejects missing, misaligned, or unavailable source chunks.
    pub fn shared(
        descriptor: ChunkDescriptor,
        structure: Structure,
        name: Arc<str>,
        storage_chunk: usize,
    ) -> Result<Self, ProviderError> {
        let (atom_start, chunk_rows, property_rows, atom_count) = {
            let chunk = structure.data().chunks.get(storage_chunk).ok_or(
                ProviderError::StorageChunkUnavailable {
                    index: storage_chunk,
                },
            )?;
            let property = structure.annotations().get(&name).ok_or_else(|| {
                ProviderError::PropertyUnavailable {
                    name: Box::from(name.as_ref()),
                }
            })?;
            (
                chunk.atoms().start,
                chunk.len(),
                property.len(),
                structure.atom_count(),
            )
        };
        if chunk_rows != descriptor.rows() {
            return Err(ProviderError::PayloadLengthMismatch {
                declared: descriptor.rows(),
                actual: chunk_rows,
            });
        }
        if property_rows != atom_count {
            return Err(ProviderError::PropertyLengthMismatch {
                name: Box::from(name.as_ref()),
                expected: atom_count,
                actual: property_rows,
            });
        }
        Ok(Self {
            descriptor,
            structure,
            name,
            atom_start,
            storage_chunk,
        })
    }

    /// Stable chunk metadata.
    #[must_use]
    pub const fn descriptor(&self) -> ChunkDescriptor {
        self.descriptor
    }

    /// Exact property name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Physical property type.
    #[must_use]
    pub fn kind(&self) -> Option<PropertyKind> {
        Some(match self.property()? {
            AtomAnnotation::Boolean(_) => PropertyKind::Boolean,
            AtomAnnotation::Integer(_) => PropertyKind::Integer,
            AtomAnnotation::Real(_) => PropertyKind::Real,
            AtomAnnotation::Symbol(_) => PropertyKind::Symbol,
        })
    }

    /// Reads one local value without allocating or decoding unrelated rows.
    ///
    /// # Errors
    ///
    /// Rejects a local row outside this chunk.
    pub fn value(&self, local: LocalRow) -> Result<PropertyValue, ProviderError> {
        let _ = self.descriptor.resolve(local)?;
        let atom =
            self.atom_start
                .checked_add(local.get())
                .ok_or(ProviderError::IdentityOverflow {
                    dataset: self.descriptor.dataset(),
                    field: "structure-local atom row",
                })?;
        let property = self
            .property()
            .ok_or_else(|| ProviderError::PropertyUnavailable {
                name: Box::from(self.name.as_ref()),
            })?;
        match property {
            AtomAnnotation::Boolean(column) => column
                .get(atom)
                .map(|(value, presence)| PropertyValue::Boolean(value, presence)),
            AtomAnnotation::Integer(column) => column
                .get(atom)
                .map(|(value, presence)| PropertyValue::Integer(value, presence)),
            AtomAnnotation::Real(column) => column
                .get(atom)
                .map(|(value, presence)| PropertyValue::Real(value, presence)),
            AtomAnnotation::Symbol(column) => column
                .get(atom)
                .map(|(value, presence)| PropertyValue::Symbol(value, presence)),
        }
        .ok_or(ProviderError::LocalRowOutOfRange {
            row: local.get(),
            rows: self.descriptor.rows(),
        })
    }

    /// Borrowed integer values when this property has integer storage.
    #[must_use]
    pub fn integers(&self) -> Option<&[i64]> {
        let Some(AtomAnnotation::Integer(column)) = self.property() else {
            return None;
        };
        slice_values(column.values(), self.atom_start, self.descriptor.rows())
    }

    /// Borrowed real values when this property has real storage.
    #[must_use]
    pub fn reals(&self) -> Option<&[f64]> {
        let Some(AtomAnnotation::Real(column)) = self.property() else {
            return None;
        };
        slice_values(column.values(), self.atom_start, self.descriptor.rows())
    }

    /// Conservative bounds retained in the source chunk metadata.
    #[must_use]
    pub fn bounds(&self) -> crate::Aabb {
        match self.structure.data().chunks.get(self.storage_chunk) {
            Some(chunk) => chunk.stats().bounds,
            None => crate::Aabb::EMPTY,
        }
    }

    fn property(&self) -> Option<&AtomAnnotation> {
        self.structure.annotations().get(&self.name)
    }
}

/// One coordinate-frame chunk retaining shared aligned storage.
#[derive(Clone, Debug)]
pub struct FrameChunk {
    descriptor: ChunkDescriptor,
    coordinates: CoordinateBlock,
    range: Range<u32>,
}

impl FrameChunk {
    /// Retains a shared coordinate block and exposes one bounded range.
    ///
    /// # Errors
    ///
    /// Rejects invalid ranges and descriptor/payload length mismatches.
    pub fn shared(
        descriptor: ChunkDescriptor,
        coordinates: CoordinateBlock,
        range: Range<u32>,
    ) -> Result<Self, ProviderError> {
        let actual =
            range
                .end
                .checked_sub(range.start)
                .ok_or(ProviderError::PayloadLengthMismatch {
                    declared: descriptor.rows(),
                    actual: 0,
                })?;
        if actual != descriptor.rows() || coordinates.range(range.clone()).is_none() {
            return Err(ProviderError::PayloadLengthMismatch {
                declared: descriptor.rows(),
                actual,
            });
        }
        Ok(Self {
            descriptor,
            coordinates,
            range,
        })
    }

    /// Stable chunk metadata.
    #[must_use]
    pub const fn descriptor(&self) -> ChunkDescriptor {
        self.descriptor
    }

    /// Coordinates borrowed directly from shared aligned storage.
    #[must_use]
    pub fn positions(&self) -> &[[f32; 3]] {
        match self.coordinates.range(self.range.clone()) {
            Some(positions) => positions,
            None => &[],
        }
    }

    /// Shared backing block for another zero-copy adapter.
    #[must_use]
    pub const fn coordinates(&self) -> &CoordinateBlock {
        &self.coordinates
    }

    /// Conservative bounds computed over only this resident chunk.
    #[must_use]
    pub fn bounds(&self) -> crate::Aabb {
        self.coordinates.bounds(self.range.clone())
    }
}

fn slice_values<T>(values: &[T], start: u32, rows: u32) -> Option<&[T]> {
    let start = usize::try_from(start).ok()?;
    let rows = usize::try_from(rows).ok()?;
    let end = start.checked_add(rows)?;
    values.get(start..end)
}

#[cfg(test)]
#[path = "payload_tests.rs"]
mod tests;
