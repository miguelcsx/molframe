use super::{parse_mol2_record, write_mol2};

const RECORD: &str = "@<TRIPOS>MOLECULE\nligand\n2 1 1 0 0\nSMALL\nUSER_CHARGES\nSYSTEM\ncomment\n@<TRIPOS>ATOM\n10 C1 0.0 1.0 2.0 C.ar 1 LIG -0.125 BACKBONE\n20 N1 1.5 1.0 2.0 N.am 1 LIG 0.125\n@<TRIPOS>BOND\n7 10 20 ar TYPE1\n@<TRIPOS>SUBSTRUCTURE\n1 LIG 1 GROUP\n";

#[test]
fn roundtrips_tripose_metadata_and_extra_sections() {
    let record = parse_mol2_record(RECORD).expect("parse MOL2");
    assert_eq!(record.atom_metadata[0].id, 10);
    assert_eq!(record.atom_metadata[0].charge, Some(-0.125));
    assert_eq!(record.bond_metadata[0].bond_type.as_ref(), "ar");
    assert_eq!(record.extra_sections[0].name.as_ref(), "SUBSTRUCTURE");
    let encoded = write_mol2(&record).expect("write MOL2");
    assert_eq!(parse_mol2_record(&encoded).expect("reparse MOL2"), record);
}

#[test]
fn complete_parser_exposes_graph_and_metadata_together() {
    let record = parse_mol2_record(RECORD).expect("MOL2 graph");
    assert_eq!(record.molecule.atoms.len(), 2);
    assert_eq!(record.molecule.bonds[0].order, 4);
}

#[test]
fn rejects_unresolved_bond_endpoint() {
    let invalid = RECORD.replace("7 10 20 ar", "7 10 99 ar");
    assert!(parse_mol2_record(&invalid).is_err());
}
