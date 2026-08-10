//! Arrow views of residue, chain and bond tables.

use crate::extension::{ExportCost, field};
use crate::stream::{ArrowStream, ArrowTableExport};
use arrow::array::{ArrayRef, Int32Array, UInt8Array, UInt32Array};
use arrow::datatypes::{DataType, Schema, SchemaRef};
use arrow::error::{ArrowError, Result};
use arrow::record_batch::RecordBatch;
use pdbiox_core::{BondOrder, BondProvenance, ChainIndex, ResidueIndex, Structure};
use std::sync::Arc;

use crate::numeric::usize_to_u32;

macro_rules! table_type {
    ($name:ident, $schema:ident, $batch:ident, $description:literal) => {
        #[doc = $description]
        #[derive(Clone, Debug)]
        pub struct $name {
            structure: Structure,
            schema: SchemaRef,
        }

        impl $name {
            /// Binds an immutable structure snapshot to this Arrow table.
            #[must_use]
            pub fn new(structure: &Structure) -> Self {
                Self {
                    structure: structure.clone(),
                    schema: Arc::new($schema()),
                }
            }

            /// Arrow schema with pdbiox extension and copy-cost metadata.
            #[must_use]
            pub fn schema(&self) -> SchemaRef {
                self.schema.clone()
            }

            /// Materialises the derived table as one record batch.
            ///
            /// # Errors
            ///
            /// Returns an Arrow error if columns violate the declared schema.
            pub fn record_batches(&self) -> Result<Vec<RecordBatch>> {
                <Self as ArrowTableExport>::record_batches(self)
            }

            /// Creates a consumable Arrow C Stream Interface value.
            ///
            /// # Errors
            ///
            /// Returns the same errors as `record_batches`.
            pub fn arrow_stream(&self) -> Result<ArrowStream> {
                <Self as ArrowTableExport>::arrow_stream(self)
            }
        }

        impl ArrowTableExport for $name {
            fn schema(&self) -> SchemaRef {
                self.schema.clone()
            }

            fn record_batches(&self) -> Result<Vec<RecordBatch>> {
                Ok(vec![$batch(&self.structure, self.schema.clone())?])
            }
        }
    };
}

table_type!(
    ResidueTable,
    residue_schema,
    residue_batch,
    "Arrow view of the residue table."
);
table_type!(
    ChainTable,
    chain_schema,
    chain_batch,
    "Arrow view of the chain table."
);
table_type!(
    BondTable,
    bond_schema,
    bond_batch,
    "Arrow view of the chemical bond table."
);

fn residue_batch(structure: &Structure, schema: SchemaRef) -> Result<RecordBatch> {
    let residues = &structure.data().topology.residues;
    let positions = 0..usize_to_u32(residues.len());
    let index = UInt32Array::from_iter_values(positions.clone());
    let chain = UInt32Array::from_iter_values(
        positions
            .clone()
            .filter_map(|position| structure.data().topology.chains.containing(position))
            .map(ChainIndex::get),
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
            .map(pdbiox_core::SymbolId::get)
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

fn chain_batch(structure: &Structure, schema: SchemaRef) -> Result<RecordBatch> {
    let chains = &structure.data().topology.chains;
    let positions = 0..usize_to_u32(chains.len());
    let index = UInt32Array::from_iter_values(positions.clone());
    let label = UInt32Array::from_iter_values(positions.clone().filter_map(|position| {
        chains
            .label_asym_id(ChainIndex::new(position))
            .map(pdbiox_core::SymbolId::get)
    }));
    let entity = UInt32Array::from_iter_values(positions.clone().filter_map(|position| {
        chains
            .entity(ChainIndex::new(position))
            .map(pdbiox_core::EntityIndex::get)
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

fn range_width(range: std::ops::Range<u32>) -> Result<u32> {
    range.end.checked_sub(range.start).ok_or_else(|| {
        ArrowError::InvalidArgumentError("topology range ends before it starts".to_owned())
    })
}

fn bond_batch(structure: &Structure, schema: SchemaRef) -> Result<RecordBatch> {
    let bonds = &structure.data().bonds;
    let index = UInt32Array::from_iter_values(0..usize_to_u32(bonds.len()));
    let atom_a = UInt32Array::from_iter_values(bonds.iter().map(|bond| bond.atom_a.get()));
    let atom_b = UInt32Array::from_iter_values(bonds.iter().map(|bond| bond.atom_b.get()));
    let order = UInt8Array::from_iter_values(bonds.iter().map(|bond| order_code(bond.order)));
    let provenance =
        UInt8Array::from_iter_values(bonds.iter().map(|bond| provenance_code(bond.provenance)));
    RecordBatch::try_new(
        schema,
        vec![
            Arc::new(index) as ArrayRef,
            Arc::new(atom_a),
            Arc::new(atom_b),
            Arc::new(order),
            Arc::new(provenance),
        ],
    )
}

fn residue_schema() -> Schema {
    Schema::new(vec![
        decoded(
            "residue_index",
            DataType::UInt32,
            "pdbiox.residue_index",
            false,
        ),
        decoded("chain_index", DataType::UInt32, "pdbiox.chain_index", false),
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
        decoded("label_comp_id", DataType::UInt32, "pdbiox.symbol_id", false),
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

fn chain_schema() -> Schema {
    Schema::new(vec![
        decoded("chain_index", DataType::UInt32, "pdbiox.chain_index", false),
        decoded("label_asym_id", DataType::UInt32, "pdbiox.symbol_id", false),
        decoded(
            "entity_index",
            DataType::UInt32,
            "pdbiox.entity_index",
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

fn bond_schema() -> Schema {
    Schema::new(vec![
        field(
            "bond_index",
            DataType::UInt32,
            false,
            None,
            ExportCost::Decode,
        ),
        decoded("atom_a", DataType::UInt32, "pdbiox.atom_index", false),
        decoded("atom_b", DataType::UInt32, "pdbiox.atom_index", false),
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

fn decoded(
    name: &str,
    data_type: DataType,
    extension: &str,
    nullable: bool,
) -> arrow::datatypes::Field {
    field(
        name,
        data_type,
        nullable,
        Some(extension),
        ExportCost::Decode,
    )
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

#[cfg(test)]
#[path = "tables_tests.rs"]
mod tests;
