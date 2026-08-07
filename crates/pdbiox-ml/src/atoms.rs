//! Arrow atom-table batches aligned to pdbiox chunk boundaries.

use crate::extension::{ExportCost, field};
use crate::owner::{f32_buffer, symbol_buffer};
use crate::stream::{ArrowStream, ArrowTableExport};
use arrow::array::{
    ArrayRef, BooleanArray, FixedSizeListArray, Float32Array, UInt8Array, UInt32Array,
};
use arrow::buffer::NullBuffer;
use arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use arrow::error::{ArrowError, Result};
use arrow::record_batch::RecordBatch;
use pdbiox_core::{AtomChunk, Structure};
use std::sync::Arc;

/// Arrow view of a structure's atom table.
#[derive(Clone, Debug)]
pub struct AtomTable {
    structure: Structure,
    schema: SchemaRef,
}

impl AtomTable {
    /// Binds the immutable structure snapshot to an Arrow table.
    #[must_use]
    pub fn new(structure: &Structure) -> Self {
        Self {
            structure: structure.clone(),
            schema: Arc::new(atom_schema()),
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
    /// # Errors
    ///
    /// Returns the same errors as [`AtomTable::record_batches`].
    pub fn arrow_stream(&self) -> Result<ArrowStream> {
        <Self as ArrowTableExport>::arrow_stream(self)
    }

    fn batch(&self, chunk: &AtomChunk) -> Result<RecordBatch> {
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
        let residue = UInt32Array::from_iter_values((0..chunk.len()).filter_map(|local| {
            chunk
                .residue(local, &self.structure.data().topology.residues)
                .map(pdbiox_core::ResidueIndex::get)
        }));
        let element = UInt8Array::from_iter_values((0..chunk.len()).filter_map(|local| {
            chunk
                .element(local)
                .map(pdbiox_core::Element::atomic_number)
        }));
        let coordinates = coordinates_array(positions, chunk, &self.structure)?;
        let names: ArrayRef = match chunk.atom_names_plain() {
            Some(names) => Arc::new(UInt32Array::new(
                symbol_buffer(names, &self.structure)?,
                None,
            )),
            None => Arc::new(UInt32Array::from_iter_values((0..chunk.len()).filter_map(
                |local| chunk.atom_name(local).map(pdbiox_core::SymbolId::get),
            ))),
        };
        let altloc = UInt32Array::from_iter_values(
            (0..chunk.len()).filter_map(|local| chunk.alt_id(local).map(pdbiox_core::AltId::get)),
        );
        let occupancy = Float32Array::new(
            f32_buffer(occupancies, &self.structure)?,
            validity(chunk.len(), |local| chunk.occupancy(local)),
        );
        let b_factor = Float32Array::new(
            f32_buffer(b_factors, &self.structure)?,
            validity(chunk.len(), |local| chunk.b_factor(local)),
        );
        let arrays: Vec<ArrayRef> = vec![
            Arc::new(index),
            Arc::new(residue),
            Arc::new(coordinates),
            Arc::new(element),
            names,
            Arc::new(altloc),
            Arc::new(occupancy),
            Arc::new(unknown_flags(chunk.len(), |local| chunk.occupancy(local))),
            Arc::new(b_factor),
            Arc::new(unknown_flags(chunk.len(), |local| chunk.b_factor(local))),
        ];
        RecordBatch::try_new(self.schema.clone(), arrays)
    }
}

impl ArrowTableExport for AtomTable {
    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }

    fn record_batches(&self) -> Result<Vec<RecordBatch>> {
        self.structure
            .data()
            .chunks
            .iter()
            .map(|chunk| self.batch(chunk))
            .collect()
    }
}

fn coordinates_array(
    positions: &[[f32; 3]],
    chunk: &AtomChunk,
    owner: &Structure,
) -> Result<FixedSizeListArray> {
    let count = positions
        .len()
        .checked_mul(3)
        .ok_or_else(|| ArrowError::MemoryError("coordinate length overflow".to_owned()))?;
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
    let flags: Vec<bool> = (0..len)
        .map(|position| value(position).is_some_and(|(_, state)| state.is_present()))
        .collect();
    (!flags.iter().all(|flag| *flag)).then(|| NullBuffer::from(flags))
}

fn unknown_flags<T>(
    len: u32,
    mut value: impl FnMut(u32) -> Option<(T, pdbiox_core::Presence)>,
) -> BooleanArray {
    (0..len)
        .map(|position| {
            value(position).is_some_and(|(_, state)| state == pdbiox_core::Presence::Unknown)
        })
        .collect()
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
