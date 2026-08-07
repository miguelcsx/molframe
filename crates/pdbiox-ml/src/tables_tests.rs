use super::*;
use pdbiox_core::io::{InputBuffer, ReadOptions};

const SOURCE: &str = "data_t\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA ALA A 1 1 2 3\nATOM 2 N N ALA A 1 4 5 6\n";

fn structure() -> Structure {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    match pdbiox_cif::read(&input, &ReadOptions::new()) {
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

fn rows(result: Result<Vec<RecordBatch>>) -> usize {
    match result {
        Ok(batches) => batches.first().map_or(0, RecordBatch::num_rows),
        Err(error) => panic!("Arrow export failed: {error}"),
    }
}
