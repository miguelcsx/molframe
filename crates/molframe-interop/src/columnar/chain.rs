//! Arrow view of the chain table.
//!
//! One record batch per fixed table batch row count; each batch materialises
//! its five columns in a single pass over the batch's chain range.

use super::extension::{ExportCost, field};
use super::table::{decoded, range_width, table_type};
use crate::numeric::usize_to_u32;
use arrow::array::{ArrayRef, UInt32Array};
use arrow::datatypes::{DataType, Schema, SchemaRef};
use arrow::error::Result;
use arrow::record_batch::RecordBatch;
use molframe_core::{ChainIndex, EntityIndex, Structure};
use std::ops::Range;
use std::sync::Arc;

table_type!(
    ChainTable,
    chain_schema,
    chain_count,
    chain_batch,
    "Arrow view of the chain table."
);

fn chain_count(structure: &Structure) -> usize {
    structure.data().topology.chains.len()
}

fn chain_batch(
    structure: &Structure,
    schema: SchemaRef,
    range: Range<usize>,
) -> Result<RecordBatch> {
    let chains = &structure.data().topology.chains;
    let positions = usize_to_u32(range.start)..usize_to_u32(range.end);
    let index = UInt32Array::from_iter_values(positions.clone());
    let label = UInt32Array::from_iter_values(positions.clone().filter_map(|position| {
        chains
            .label_asym_id(ChainIndex::new(position))
            .map(molframe_core::SymbolId::get)
    }));
    let entity = UInt32Array::from_iter_values(positions.clone().filter_map(|position| {
        chains
            .entity(ChainIndex::new(position))
            .map(EntityIndex::get)
    }));
    let first_residue = UInt32Array::from_iter_values(
        positions
            .clone()
            .filter_map(|position| chains.residues(ChainIndex::new(position)))
            .map(|range| range.start),
    );
    let residue_count = positions
        .filter_map(|position| chains.residues(ChainIndex::new(position)))
        .map(range_width)
        .collect::<Result<Vec<_>>>()?;
    let residue_count = UInt32Array::from(residue_count);
    RecordBatch::try_new(
        schema,
        vec![
            Arc::new(index) as ArrayRef,
            Arc::new(label),
            Arc::new(entity),
            Arc::new(first_residue),
            Arc::new(residue_count),
        ],
    )
}

fn chain_schema() -> Schema {
    Schema::new(vec![
        decoded(
            "chain_index",
            DataType::UInt32,
            "molframe.chain_index",
            false,
        ),
        decoded(
            "label_asym_id",
            DataType::UInt32,
            "molframe.symbol_id",
            false,
        ),
        decoded(
            "entity_index",
            DataType::UInt32,
            "molframe.entity_index",
            false,
        ),
        field(
            "first_residue",
            DataType::UInt32,
            false,
            None,
            ExportCost::Decode,
        ),
        field(
            "residue_count",
            DataType::UInt32,
            false,
            None,
            ExportCost::Decode,
        ),
    ])
}
