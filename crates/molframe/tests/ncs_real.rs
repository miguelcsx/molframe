//! NCS operators of deposited entries bind to the chains they relate, and
//! survive a write and read in every text and binary format.

use molframe::ReadOptions;
use molframe::crystal::{NcsBindOptions, NcsCode, NcsExt, NcsSet};

const GROEL_GZ: &[u8] = include_bytes!("../../molframe-bench/data/1aon.cif.gz");
const FOUR_HHB_PDB: &[u8] = include_bytes!("../../molframe-bench/data/4hhb.pdb");

fn read(bytes: &[u8], name: &str) -> molframe::Structure {
    match molframe::read_bytes(bytes.to_vec(), Some(name), &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("{name} must read: {findings:?}"),
    }
}

/// One operator flattened to `(id, code, rotation, translation)`.
type Operator = (String, NcsCode, [[f64; 3]; 3], [f64; 3]);

fn operators(set: &NcsSet) -> Vec<Operator> {
    set.operators()
        .map(|operator| {
            (
                operator.id.to_string(),
                operator.code,
                operator.transform.matrix,
                operator.transform.translation,
            )
        })
        .collect()
}

fn assert_same(before: &NcsSet, after: &NcsSet, tolerance: f64) {
    let (before, after) = (operators(before), operators(after));
    assert_eq!(before.len(), after.len());
    for (left, right) in before.iter().zip(&after) {
        assert_eq!(left.1, right.1);
        for row in 0..3 {
            assert!((left.3[row] - right.3[row]).abs() <= tolerance);
            for column in 0..3 {
                assert!((left.2[row][column] - right.2[row][column]).abs() <= tolerance);
            }
        }
    }
}

#[test]
fn groel_operators_bind_to_chains_that_superpose() {
    let structure = read(GROEL_GZ, "1aon.cif.gz");
    let pairs = structure
        .ncs_chain_pairs(&NcsBindOptions::default())
        .expect("binding stays inside its limits");
    assert!(!pairs.is_empty(), "the deposited operators relate chains");
    for pair in &pairs {
        assert!(pair.rmsd.is_some_and(|value| value <= 2.0));
        assert_ne!(pair.source, pair.target);
    }
}

#[test]
fn groel_operators_survive_cif_bcif_and_pdb_writes() {
    let structure = read(GROEL_GZ, "1aon.cif.gz");
    let original = structure.ncs_set().expect("NCS attached on read");
    assert_eq!(original.len(), 18);

    let cif = molframe::write_mmcif(&structure).expect("mmCIF write");
    let again = read(cif.as_bytes(), "again.cif");
    assert_same(original, again.ncs_set().expect("NCS in mmCIF"), 1e-5);

    let bcif = molframe::write_bcif(&structure).expect("BinaryCIF write");
    let again = read(&bcif, "again.bcif");
    assert_same(original, again.ncs_set().expect("NCS in BinaryCIF"), 1e-5);

    let options = molframe::formats::pdb::PdbOptions::new()
        .namespace(molframe::formats::pdb::PdbIdentifierNamespace::Auth);
    let pdb = molframe::write_pdb(&structure, &options).expect("PDB write");
    let again = read(pdb.as_bytes(), "again.pdb");
    assert_same(original, again.ncs_set().expect("NCS in PDB"), 1e-5);
}

#[test]
fn a_pdb_entry_with_the_identity_operator_binds_nothing_across_chains() {
    let structure = read(FOUR_HHB_PDB, "4hhb.pdb");
    let pairs = structure
        .ncs_chain_pairs(&NcsBindOptions::default())
        .expect("binding stays inside its limits");
    assert!(pairs.iter().all(|pair| pair.source != pair.target));
    let expanded = structure.ncs_expanded().expect("no generate operator");
    assert_eq!(expanded.atom_count(), structure.atom_count());
}
