use super::table::TABLE_BATCH_ROWS;
use super::*;
use arrow::error::Result;
use arrow::record_batch::RecordBatch;
use molframe_core::Structure;
use molframe_core::io::{InputBuffer, ReadOptions};

const SOURCE: &str = "data_t\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA ALA A 1 1 2 3\nATOM 2 N N ALA A 1 4 5 6\n";

fn structure() -> Structure {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

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
        });
    }
    data.bonds = bonds.finish();
    let table = BondTable::new(&Structure::new(data));
    let batches = table.record_batches().expect("bounded bond batches");
    assert_eq!(batches.len(), 2);
    assert_eq!(batches[0].num_rows(), TABLE_BATCH_ROWS);
    assert_eq!(batches[1].num_rows(), 1);
}

fn rows(result: Result<Vec<RecordBatch>>) -> usize {
    match result {
        Ok(batches) => batches.iter().map(RecordBatch::num_rows).sum(),
        Err(error) => panic!("Arrow export failed: {error}"),
    }
}
