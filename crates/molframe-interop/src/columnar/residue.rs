//! Arrow view of the residue table.
//!
//! One record batch per fixed table batch row count; each batch materialises
//! its seven columns in a single pass over the batch's residue range.

use super::extension::{ExportCost, field};
use super::table::{decoded, range_width, table_type};
use crate::numeric::usize_to_u32;
use arrow::array::{ArrayRef, Int32Array, UInt32Array};
use arrow::datatypes::{DataType, Schema, SchemaRef};
use arrow::error::Result;
use arrow::record_batch::RecordBatch;
use molframe_core::{ResidueIndex, Structure, SymbolId};
use std::ops::Range;
use std::sync::Arc;

table_type!(
    ResidueTable,
    residue_schema,
    residue_count,
    residue_batch,
    "Arrow view of the residue table."
);

fn residue_count(structure: &Structure) -> usize {
    structure.data().topology.residues.len()
}

fn residue_batch(
    structure: &Structure,
    schema: SchemaRef,
    range: Range<usize>,
) -> Result<RecordBatch> {
    let residues = &structure.data().topology.residues;
    let positions = usize_to_u32(range.start)..usize_to_u32(range.end);
    let index = UInt32Array::from_iter_values(positions.clone());
    let chain = UInt32Array::from_iter_values(
        positions
            .clone()
            .filter_map(|position| structure.data().topology.chains.containing(position))
            .map(molframe_core::ChainIndex::get),
    );
    let first_atom = UInt32Array::from_iter_values(
        positions
            .clone()
            .filter_map(|position| residues.atoms(ResidueIndex::new(position)))
            .map(|range| range.start),
    );
    let atom_count = positions
        .clone()
        .filter_map(|position| residues.atoms(ResidueIndex::new(position)))
        .map(range_width)
        .collect::<Result<Vec<_>>>()?;
    let atom_count = UInt32Array::from(atom_count);
    let component = UInt32Array::from_iter_values(positions.clone().filter_map(|position| {
        residues
            .label_comp_id(ResidueIndex::new(position))
            .map(SymbolId::get)
    }));
    let label_sequence = positions
        .clone()
        .map(|position| residues.label_seq_id(ResidueIndex::new(position)))
        .collect::<Int32Array>();
    let auth_sequence = positions
        .map(|position| residues.auth_seq_id(ResidueIndex::new(position)))
        .collect::<Int32Array>();
    RecordBatch::try_new(
        schema,
        vec![
            Arc::new(index) as ArrayRef,
            Arc::new(chain),
            Arc::new(first_atom),
            Arc::new(atom_count),
            Arc::new(component),
            Arc::new(label_sequence),
            Arc::new(auth_sequence),
        ],
    )
}

fn residue_schema() -> Schema {
    Schema::new(vec![
        decoded(
            "residue_index",
            DataType::UInt32,
            "molframe.residue_index",
            false,
        ),
        decoded(
            "chain_index",
            DataType::UInt32,
            "molframe.chain_index",
            false,
        ),
        field(
            "first_atom",
            DataType::UInt32,
            false,
            None,
            ExportCost::Decode,
        ),
        field(
            "atom_count",
            DataType::UInt32,
            false,
            None,
            ExportCost::Decode,
        ),
        decoded(
            "label_comp_id",
            DataType::UInt32,
            "molframe.symbol_id",
            false,
        ),
        field(
            "label_seq_id",
            DataType::Int32,
            true,
            None,
            ExportCost::Decode,
        ),
        field(
            "auth_seq_id",
            DataType::Int32,
            true,
            None,
            ExportCost::Decode,
        ),
    ])
}
