//! Arrow view of the chemical bond table.
//!
//! One record batch per fixed table batch row count; each batch materialises
//! its five columns in a single pass over the batch's bond range.

use super::extension::{ExportCost, field};
use super::table::{decoded, table_type};
use crate::numeric::usize_to_u32;
use arrow::array::{ArrayRef, UInt8Array, UInt32Array};
use arrow::datatypes::{DataType, Schema, SchemaRef};
use arrow::error::{ArrowError, Result};
use arrow::record_batch::RecordBatch;
use molframe_core::{BondIndex, BondOrder, BondProvenance, Structure};
use std::ops::Range;
use std::sync::Arc;

table_type!(
    BondTable,
    bond_schema,
    bond_count,
    bond_batch,
    "Arrow view of the chemical bond table."
);

fn bond_count(structure: &Structure) -> usize {
    structure.data().bonds.len()
}

fn bond_batch(
    structure: &Structure,
    schema: SchemaRef,
    range: Range<usize>,
) -> Result<RecordBatch> {
    let bonds = &structure.data().bonds;
    let start = usize_to_u32(range.start);
    let end = usize_to_u32(range.end);
    let mut index = Vec::with_capacity(range.len());
    let mut atom_a = Vec::with_capacity(range.len());
    let mut atom_b = Vec::with_capacity(range.len());
    let mut order = Vec::with_capacity(range.len());
    let mut provenance = Vec::with_capacity(range.len());
    for raw in start..end {
        let Some(bond) = bonds.get(BondIndex::new(raw)) else {
            return Err(ArrowError::InvalidArgumentError(
                "bond row is absent".to_owned(),
            ));
        };
        index.push(raw);
        atom_a.push(bond.atom_a.get());
        atom_b.push(bond.atom_b.get());
        order.push(order_code(bond.order));
        provenance.push(provenance_code(bond.provenance));
    }
    RecordBatch::try_new(
        schema,
        vec![
            Arc::new(UInt32Array::from(index)) as ArrayRef,
            Arc::new(UInt32Array::from(atom_a)),
            Arc::new(UInt32Array::from(atom_b)),
            Arc::new(UInt8Array::from(order)),
            Arc::new(UInt8Array::from(provenance)),
        ],
    )
}

fn bond_schema() -> Schema {
    Schema::new(vec![
        field(
            "bond_index",
            DataType::UInt32,
            false,
            None,
            ExportCost::Decode,
        ),
        decoded("atom_a", DataType::UInt32, "molframe.atom_index", false),
        decoded("atom_b", DataType::UInt32, "molframe.atom_index", false),
        field("order", DataType::UInt8, false, None, ExportCost::Decode),
        field(
            "provenance",
            DataType::UInt8,
            false,
            None,
            ExportCost::Decode,
        ),
    ])
}

const fn order_code(order: BondOrder) -> u8 {
    match order {
        BondOrder::Unknown => 0,
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Aromatic => 5,
        BondOrder::Polymeric => 6,
    }
}

const fn provenance_code(provenance: BondProvenance) -> u8 {
    match provenance {
        BondProvenance::File => 1,
        BondProvenance::ChemicalComponentDictionary => 2,
        BondProvenance::InferredDistance => 3,
        BondProvenance::User => 4,
    }
}
