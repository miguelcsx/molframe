use super::table::TABLE_BATCH_ROWS;
use super::*;
use crate::test_support::two_atoms as structure;
use arrow::error::Result;
use arrow::record_batch::RecordBatch;
use molframe_core::Structure;
use molframe_core::optional::OptionalSymbol;

#[test]
fn every_normalised_table_exports_a_c_stream() {
    let structure = structure();
    let residues = ResidueTable::new(&structure);
    let chains = ChainTable::new(&structure);
    let bonds = BondTable::new(&structure);
    assert_eq!(rows(residues.record_batches()), 1);
    assert_eq!(rows(chains.record_batches()), 1);
    assert_eq!(rows(bonds.record_batches()), 0);
    assert!(residues.arrow_stream().is_ok());
    assert!(chains.arrow_stream().is_ok());
    assert!(bonds.arrow_stream().is_ok());
}

#[test]
fn topology_exports_are_split_into_bounded_batches() {
    let mut data = molframe_core::StructureData::empty();
    let mut bonds = molframe_core::BondTableBuilder::new();
    let batch_rows = u32::try_from(TABLE_BATCH_ROWS).expect("batch rows fit u32");
    for atom in 0..=batch_rows {
        bonds.push(molframe_core::BondRecord {
            atom_a: molframe_core::AtomIndex::new(atom),
            atom_b: molframe_core::AtomIndex::new(atom + 1),
            order: molframe_core::BondOrder::Single,
            provenance: molframe_core::BondProvenance::User,
            kind: OptionalSymbol::NONE,
        });
    }
    data.bonds = bonds.finish();
    let table = BondTable::new(&Structure::new(data));
    let batches = table.record_batches().expect("bounded bond batches");
    assert_eq!(batches.len(), 2);
    assert_eq!(batches[0].num_rows(), TABLE_BATCH_ROWS);
    assert_eq!(batches[1].num_rows(), 1);
}

#[test]
fn a_row_window_over_derived_batches_slices_only_the_overlapping_ones() {
    let mut data = molframe_core::StructureData::empty();
    let mut bonds = molframe_core::BondTableBuilder::new();
    let batch_rows = u32::try_from(TABLE_BATCH_ROWS).expect("batch rows fit u32");
    for atom in 0..=batch_rows {
        bonds.push(molframe_core::BondRecord {
            atom_a: molframe_core::AtomIndex::new(atom),
            atom_b: molframe_core::AtomIndex::new(atom + 1),
            order: molframe_core::BondOrder::Single,
            provenance: molframe_core::BondProvenance::User,
            kind: OptionalSymbol::NONE,
        });
    }
    data.bonds = bonds.finish();
    let table = BondTable::new(&Structure::new(data));
    let stream = table
        .arrow_stream_rows(TABLE_BATCH_ROWS - 3..TABLE_BATCH_ROWS + 1)
        .expect("a window inside the table");
    let reader = arrow::ffi_stream::ArrowArrayStreamReader::try_new(stream.into_ffi())
        .expect("the stream imports");
    let mut firsts = Vec::new();
    for batch in reader {
        let batch = batch.expect("the window batch converts");
        let column = batch
            .column_by_name("bond_index")
            .and_then(|array| array.as_any().downcast_ref::<arrow::array::UInt32Array>())
            .expect("index column");
        firsts.extend(column.values().iter().copied());
    }
    let expected: Vec<u32> = (batch_rows - 3..=batch_rows).collect();
    assert_eq!(firsts, expected);
    assert!(table.arrow_stream_rows(0..TABLE_BATCH_ROWS + 2).is_err());
}

fn rows(result: Result<Vec<RecordBatch>>) -> usize {
    match result {
        Ok(batches) => batches.iter().map(RecordBatch::num_rows).sum(),
        Err(error) => panic!("Arrow export failed: {error}"),
    }
}
