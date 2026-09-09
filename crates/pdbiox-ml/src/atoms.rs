//! Arrow atom-table batches aligned to pdbiox chunk boundaries.

use crate::extension::{ExportCost, field};
use crate::owner::{SnapshotOwner, f32_buffer, symbol_buffer};
use crate::stream::{ArrowStream, ArrowTableExport};
use arrow::array::{
    ArrayRef, BooleanArray, BooleanBuilder, FixedSizeListArray, Float32Array, UInt8Array,
    UInt32Array, builder::NullBufferBuilder,
};
use arrow::buffer::NullBuffer;
use arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use arrow::error::{ArrowError, Result};
use arrow::record_batch::RecordBatch;
use pdbiox_core::topology::ResidueTable;
use pdbiox_core::{AtomChunk, ResidueIndex, Structure};
use std::ops::Range;
use std::sync::Arc;

/// Arrow view of a structure's atom table.
#[derive(Clone, Debug)]
pub struct AtomTable {
    structure: Structure,
    schema: SchemaRef,
    owner: SnapshotOwner,
}

impl AtomTable {
    /// Binds the immutable structure snapshot to an Arrow table.
    #[must_use]
    pub fn new(structure: &Structure) -> Self {
        Self {
            structure: structure.clone(),
            schema: Arc::new(atom_schema()),
            owner: SnapshotOwner::new(structure),
        }
    }

    /// Arrow schema with pdbiox extension and copy-cost metadata.
    #[must_use]
    pub fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }

    /// One record batch per internal atom chunk.
    ///
    /// # Errors
    ///
    /// Returns an Arrow error if a chunk violates the declared schema or a
    /// backing allocation cannot be represented.
    pub fn record_batches(&self) -> Result<Vec<RecordBatch>> {
        <Self as ArrowTableExport>::record_batches(self)
    }

    /// Creates a consumable Arrow C Stream Interface value.
    ///
    /// Construction is constant-space: chunks are converted only when the
    /// consumer pulls them from the stream. A conversion error is reported by
    /// that pull rather than by this method.
    ///
    /// # Errors
    ///
    /// Reserved for failures that occur while constructing the stream adapter.
    pub fn arrow_stream(&self) -> Result<ArrowStream> {
        <Self as ArrowTableExport>::arrow_stream(self)
    }

    fn chunk_batch(&self, chunk: &AtomChunk) -> Result<RecordBatch> {
        let range = chunk.atoms();
        let positions = self
            .structure
            .positions()
            .get(range.start as usize..range.end as usize)
            .ok_or_else(|| {
                ArrowError::InvalidArgumentError("chunk coordinates absent".to_owned())
            })?;
        let occupancies = chunk.occupancies_plain().ok_or_else(|| {
            ArrowError::InvalidArgumentError("occupancies require decode".to_owned())
        })?;
        let b_factors = chunk.b_factors_plain().ok_or_else(|| {
            ArrowError::InvalidArgumentError("temperature factors require decode".to_owned())
        })?;
        let index = UInt32Array::from_iter_values(range.clone());
        let residue = residue_indices(range.clone(), &self.structure.data().topology.residues)?;
        let element = UInt8Array::from_iter_values((0..chunk.len()).filter_map(|local| {
            chunk
                .element(local)
                .map(pdbiox_core::Element::atomic_number)
        }));
        let coordinates = coordinates_array(positions, chunk, &self.owner)?;
        let names: ArrayRef = match chunk.atom_names_plain() {
            Some(names) => Arc::new(UInt32Array::new(symbol_buffer(names, &self.owner)?, None)),
            None => Arc::new(UInt32Array::from_iter_values((0..chunk.len()).filter_map(
                |local| chunk.atom_name(local).map(pdbiox_core::SymbolId::get),
            ))),
        };
        let altloc = UInt32Array::from_iter_values(
            (0..chunk.len()).filter_map(|local| chunk.alt_id(local).map(pdbiox_core::AltId::get)),
        );
        let (occupancy_validity, occupancy_unknown) =
            presence_arrays(chunk.len(), |local| chunk.occupancy(local));
        let occupancy =
            Float32Array::new(f32_buffer(occupancies, &self.owner)?, occupancy_validity);
        let (b_factor_validity, b_factor_unknown) =
            presence_arrays(chunk.len(), |local| chunk.b_factor(local));
        let b_factor = Float32Array::new(f32_buffer(b_factors, &self.owner)?, b_factor_validity);
        let arrays: Vec<ArrayRef> = vec![
            Arc::new(index),
            Arc::new(residue),
            Arc::new(coordinates),
            Arc::new(element),
            names,
            Arc::new(altloc),
            Arc::new(occupancy),
            Arc::new(occupancy_unknown),
            Arc::new(b_factor),
            Arc::new(b_factor_unknown),
        ];
        RecordBatch::try_new(self.schema.clone(), arrays)
    }
}

/// Expands ordered residue ranges in `O(atoms + overlapping residues)` time.
fn residue_indices(range: Range<u32>, residues: &ResidueTable) -> Result<UInt32Array> {
    let atom_count = range.end.checked_sub(range.start).ok_or_else(|| {
        ArrowError::InvalidArgumentError("atom chunk range is reversed".to_owned())
    })?;
    let expected = usize::try_from(atom_count)
        .map_err(|_| ArrowError::MemoryError("atom chunk length exceeds usize".to_owned()))?;
    if expected == 0 {
        return Ok(UInt32Array::from(Vec::<u32>::new()));
    }
    let last_atom = range
        .end
        .checked_sub(1)
        .ok_or_else(|| ArrowError::InvalidArgumentError("empty atom chunk".to_owned()))?;
    let first_residue = residues.containing(range.start).ok_or_else(|| {
        ArrowError::InvalidArgumentError("first atom has no containing residue".to_owned())
    })?;
    let last_residue = residues.containing(last_atom).ok_or_else(|| {
        ArrowError::InvalidArgumentError("last atom has no containing residue".to_owned())
    })?;
    let mut values = Vec::with_capacity(expected);
    for raw in first_residue.get()..=last_residue.get() {
        let atoms = residues.atoms(ResidueIndex::new(raw)).ok_or_else(|| {
            ArrowError::InvalidArgumentError("residue atom range absent".to_owned())
        })?;
        let overlap_start = atoms.start.max(range.start);
        let overlap_end = atoms.end.min(range.end);
        let overlap = usize::try_from(overlap_end.saturating_sub(overlap_start))
            .map_err(|_| ArrowError::MemoryError("residue overlap exceeds usize".to_owned()))?;
        let next = values
            .len()
            .checked_add(overlap)
            .ok_or_else(|| ArrowError::MemoryError("residue index length overflow".to_owned()))?;
        if next > expected {
            return Err(ArrowError::InvalidArgumentError(
                "residue ranges overlap within atom chunk".to_owned(),
            ));
        }
        values.resize(next, raw);
    }
    if values.len() != expected {
        return Err(ArrowError::InvalidArgumentError(
            "residue ranges do not cover atom chunk".to_owned(),
        ));
    }
    Ok(UInt32Array::from(values))
}

impl ArrowTableExport for AtomTable {
    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }

    fn batch_count(&self) -> usize {
        self.structure.data().chunks.len()
    }

    fn batch(&self, index: usize) -> Result<RecordBatch> {
        let chunk = self.structure.data().chunks.get(index).ok_or_else(|| {
            ArrowError::InvalidArgumentError("atom chunk index is out of bounds".to_owned())
        })?;
        self.chunk_batch(chunk)
    }
}

fn coordinates_array(
    positions: &[[f32; 3]],
    chunk: &AtomChunk,
    owner: &SnapshotOwner,
) -> Result<FixedSizeListArray> {
    let count = positions
        .len()
        .checked_mul(3)
        .ok_or_else(|| ArrowError::MemoryError("coordinate length overflow".to_owned()))?;
    // SAFETY: `[f32; 3]` is three contiguous, equally aligned `f32` values with
    // no padding. `count` is the checked scalar length of the live input slice.
    let values = unsafe { std::slice::from_raw_parts(positions.as_ptr().cast::<f32>(), count) };
    let values: ArrayRef = Arc::new(Float32Array::new(f32_buffer(values, owner)?, None));
    FixedSizeListArray::try_new(
        Arc::new(Field::new("item", DataType::Float32, false)),
        3,
        values,
        validity(chunk.len(), |local| {
            chunk
                .has_position(local)
                .then_some(((), pdbiox_core::Presence::Present))
        }),
    )
}

fn validity<T>(
    len: u32,
    mut value: impl FnMut(u32) -> Option<(T, pdbiox_core::Presence)>,
) -> Option<NullBuffer> {
    let mut validity = NullBufferBuilder::new(len as usize);
    for position in 0..len {
        let present = value(position).is_some_and(|(_, state)| state.is_present());
        validity.append(present);
    }
    validity.finish()
}

fn presence_arrays<T>(
    len: u32,
    mut value: impl FnMut(u32) -> Option<(T, pdbiox_core::Presence)>,
) -> (Option<NullBuffer>, BooleanArray) {
    let mut validity = NullBufferBuilder::new(len as usize);
    let mut unknown = BooleanBuilder::with_capacity(len as usize);
    for position in 0..len {
        let state = value(position).map(|(_, state)| state);
        let present = state.is_some_and(pdbiox_core::Presence::is_present);
        validity.append(present);
        unknown.append_value(state == Some(pdbiox_core::Presence::Unknown));
    }
    (validity.finish(), unknown.finish())
}

fn atom_schema() -> Schema {
    let coordinates =
        DataType::FixedSizeList(Arc::new(Field::new("item", DataType::Float32, false)), 3);
    Schema::new(vec![
        field(
            "atom_index",
            DataType::UInt32,
            false,
            Some("pdbiox.atom_index"),
            ExportCost::Decode,
        ),
        field(
            "residue_index",
            DataType::UInt32,
            false,
            Some("pdbiox.residue_index"),
            ExportCost::Decode,
        ),
        field(
            "coordinates",
            coordinates,
            true,
            Some("pdbiox.coordinates3f"),
            ExportCost::ZeroCopy,
        ),
        field(
            "element",
            DataType::UInt8,
            false,
            Some("pdbiox.element"),
            ExportCost::Decode,
        ),
        field(
            "atom_name",
            DataType::UInt32,
            false,
            Some("pdbiox.symbol_id"),
            ExportCost::Decode,
        ),
        field(
            "altloc",
            DataType::UInt32,
            false,
            Some("pdbiox.altloc"),
            ExportCost::Decode,
        ),
        field(
            "occupancy",
            DataType::Float32,
            true,
            None,
            ExportCost::ZeroCopy,
        ),
        field(
            "occupancy_unknown",
            DataType::Boolean,
            false,
            Some("pdbiox.validity"),
            ExportCost::Decode,
        ),
        field(
            "b_factor",
            DataType::Float32,
            true,
            None,
            ExportCost::ZeroCopy,
        ),
        field(
            "b_factor_unknown",
            DataType::Boolean,
            false,
            Some("pdbiox.validity"),
            ExportCost::Decode,
        ),
    ])
}

#[cfg(test)]
#[path = "atoms_tests.rs"]
mod tests;
