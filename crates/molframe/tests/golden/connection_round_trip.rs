//! Explicit bonds to ligands survive a write and a re-read.
//!
//! 4HHB's four haem irons are bonded to histidines through `struct_conn`.
//! Haem is not a polymer, so it has no label sequence number; the written
//! connection must carry the author sequence the reader identifies it by, or
//! the re-read structure is rejected. This is the path a structure takes to a
//! browser, so a failure here is a viewer that cannot open the molecule.

use molframe::formats::cif::CifWriteOptions;

fn options() -> CifWriteOptions {
    CifWriteOptions::new()
        .with_block_id("molframe")
        .with_generated_connection_ids()
        .with_connection_type_id("covale")
}

/// RCSB's own `BinaryCIF` of 4HHB, the benchmark fixture.
const FOUR_HHB: &[u8] = include_bytes!("../../../molframe-bench/data/4hhb.bcif");

fn bonds(structure: &molframe::Structure) -> Vec<(u32, u32)> {
    structure
        .engine()
        .data()
        .bonds
        .iter()
        .map(|bond| (bond.atom_a.get(), bond.atom_b.get()))
        .collect()
}

fn read(bytes: Vec<u8>, name: &str) -> molframe::Structure {
    match molframe::read_bytes(bytes, Some(name), &molframe::ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("{name} must read: {findings:?}"),
    }
}

#[test]
fn bonds_to_a_ligand_survive_a_written_mmcif() {
    let original = read(FOUR_HHB.to_vec(), "4hhb.bcif");
    assert!(!bonds(&original).is_empty(), "4HHB carries explicit bonds");
    let text = match molframe::write_mmcif_with_options(&original, &options()) {
        Ok(text) => text,
        Err(findings) => panic!("4HHB must write: {findings:?}"),
    };
    let copy = read(text.into_bytes(), "copy.cif");
    assert_eq!(bonds(&copy), bonds(&original));
    assert_eq!(copy.coordinates().len(), original.coordinates().len());
}

#[test]
fn bonds_to_a_ligand_survive_a_written_binary_cif() {
    let original = read(FOUR_HHB.to_vec(), "4hhb.bcif");
    let binary = match molframe::write_bcif_with_options(&original, &options()) {
        Ok(binary) => binary,
        Err(findings) => panic!("4HHB must write: {findings:?}"),
    };
    let copy = read(binary, "copy.bcif");
    assert_eq!(bonds(&copy), bonds(&original));
    assert_eq!(copy.coordinates().len(), original.coordinates().len());
}
