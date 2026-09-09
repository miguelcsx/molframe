//! Shared bond payloads whose atom endpoints remain globally addressable.
//!
//! A chunk retains one immutable structure snapshot and a bond-table range.
//! Construction and random access are O(1); no coordinate or topology column
//! is copied. The consumer converts endpoints to compact local indices only
//! after locating the corresponding resident atom pages.

use super::{ChunkDescriptor, DatasetId, LocalRow, LogicalRow, ProviderError};
use crate::{BondIndex, BondOrder, BondProvenance, BondRecord, Structure};
use std::ops::Range;

/// Stable atom identity independent of chunk residency.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct AtomEndpoint {
    dataset: DatasetId,
    row: LogicalRow,
}

impl AtomEndpoint {
    /// Creates an endpoint in an atom dataset's global logical row space.
    #[must_use]
    pub const fn new(dataset: DatasetId, row: LogicalRow) -> Self {
        Self { dataset, row }
    }

    /// Atom dataset containing the endpoint.
    #[must_use]
    pub const fn dataset(self) -> DatasetId {
        self.dataset
    }

    /// Stable atom row, never a chunk-local index.
    #[must_use]
    pub const fn row(self) -> LogicalRow {
        self.row
    }
}

/// One bond projected from native storage with global atom endpoints.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BondChunkRecord {
    /// First global atom endpoint.
    pub atom_a: AtomEndpoint,
    /// Second global atom endpoint.
    pub atom_b: AtomEndpoint,
    /// Chemical order retained from the source.
    pub order: BondOrder,
    /// Provenance retained from the source.
    pub provenance: BondProvenance,
}

/// One independently requestable range of an immutable bond table.
#[derive(Clone, Debug)]
pub struct BondChunk {
    descriptor: ChunkDescriptor,
    structure: Structure,
    source_range: Range<u32>,
    atom_dataset: DatasetId,
    atom_logical_start: LogicalRow,
}

impl BondChunk {
    /// Retains an existing bond-table range without materialising edge columns.
    ///
    /// # Errors
    ///
    /// Rejects unavailable topology, mismatched lengths, ranges outside the
    /// source table, or atom identities that overflow the global row space.
    pub fn shared(
        descriptor: ChunkDescriptor,
        structure: Structure,
        source_range: Range<u32>,
        atom_dataset: DatasetId,
        atom_logical_start: LogicalRow,
    ) -> Result<Self, ProviderError> {
        if !structure.data().bonds.is_available() {
            return Err(ProviderError::BondTopologyUnavailable);
        }
        let actual = source_range.end.checked_sub(source_range.start).ok_or(
            ProviderError::BondRangeUnavailable {
                first: u64::from(source_range.start),
                end: u64::from(source_range.end),
                rows: structure.data().bonds.len() as u64,
            },
        )?;
        if actual != descriptor.rows() {
            return Err(ProviderError::PayloadLengthMismatch {
                declared: descriptor.rows(),
                actual,
            });
        }
        let table_rows = u64::try_from(structure.data().bonds.len()).map_err(|_| {
            ProviderError::IdentityOverflow {
                dataset: descriptor.dataset(),
                field: "bond table rows",
            }
        })?;
        if u64::from(source_range.end) > table_rows {
            return Err(ProviderError::BondRangeUnavailable {
                first: u64::from(source_range.start),
                end: u64::from(source_range.end),
                rows: table_rows,
            });
        }
        validate_atom_span(
            descriptor.dataset(),
            atom_logical_start,
            structure.atom_count(),
        )?;
        Ok(Self {
            descriptor,
            structure,
            source_range,
            atom_dataset,
            atom_logical_start,
        })
    }

    /// Stable bond-chunk metadata.
    #[must_use]
    pub const fn descriptor(&self) -> ChunkDescriptor {
        self.descriptor
    }

    /// Atom dataset referenced by both endpoint columns.
    #[must_use]
    pub const fn atom_dataset(&self) -> DatasetId {
        self.atom_dataset
    }

    /// Retained immutable source; its coordinate columns are never copied.
    #[must_use]
    pub const fn structure(&self) -> &Structure {
        &self.structure
    }

    /// Projects one local bond row to global atom identities in O(1).
    ///
    /// # Errors
    ///
    /// Returns a typed error for a local row outside the chunk or identity
    /// arithmetic overflow.
    pub fn record(&self, local: LocalRow) -> Result<BondChunkRecord, ProviderError> {
        self.descriptor.resolve(local)?;
        let source = self.source_range.start.checked_add(local.get()).ok_or(
            ProviderError::IdentityOverflow {
                dataset: self.descriptor.dataset(),
                field: "bond source row",
            },
        )?;
        let record = self
            .structure
            .data()
            .bonds
            .get(BondIndex::new(source))
            .ok_or(ProviderError::BondRangeUnavailable {
                first: u64::from(source),
                end: u64::from(source) + 1,
                rows: self.structure.data().bonds.len() as u64,
            })?;
        self.global_record(record)
    }

    /// Minimum bytes represented by the retained native edge columns.
    ///
    /// # Errors
    ///
    /// Returns a typed error if byte accounting exceeds `u64`.
    pub fn minimum_host_bytes(&self) -> Result<u64, ProviderError> {
        u64::from(self.descriptor.rows())
            .checked_mul(std::mem::size_of::<BondRecord>() as u64)
            .ok_or(ProviderError::IdentityOverflow {
                dataset: self.descriptor.dataset(),
                field: "bond payload bytes",
            })
    }

    fn global_record(&self, record: BondRecord) -> Result<BondChunkRecord, ProviderError> {
        Ok(BondChunkRecord {
            atom_a: self.endpoint(record.atom_a.get())?,
            atom_b: self.endpoint(record.atom_b.get())?,
            order: record.order,
            provenance: record.provenance,
        })
    }

    fn endpoint(&self, atom: u32) -> Result<AtomEndpoint, ProviderError> {
        let row = self
            .atom_logical_start
            .get()
            .checked_add(u64::from(atom))
            .ok_or(ProviderError::IdentityOverflow {
                dataset: self.atom_dataset,
                field: "atom endpoint",
            })?;
        Ok(AtomEndpoint::new(self.atom_dataset, LogicalRow::new(row)))
    }
}

fn validate_atom_span(
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

#[cfg(test)]
#[path = "bond_payload_tests.rs"]
mod tests;
