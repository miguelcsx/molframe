use super::*;
use molframe_core::diagnostic::Code;
use molframe_core::io::InputBuffer;

const HEAD: &str = "data_s\n\
loop_\n_chem_comp.id\n_chem_comp.name\n_chem_comp.type\n\
LIG 'a ligand' NON-POLYMER\nOTH 'another' NON-POLYMER\n\
loop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n\
_chem_comp_atom.pdbx_aromatic_flag\n";
const BOND_HEAD: &str = "loop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n\
_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n_chem_comp_bond.pdbx_aromatic_flag\n";

fn document(atoms: &str, bonds: &str) -> Document {
    let text = format!("{HEAD}{atoms}{BOND_HEAD}{bonds}");
    let input = InputBuffer::from_bytes(text.into_bytes());
    molframe_cif::parse(&input).expect("parses").0
}

fn provider(document: &Document) -> Result<EmbeddedProvider, Vec<Diagnostic>> {
    EmbeddedProvider::from_document(document, DictionaryVersion::new("t"))
}

fn first_code(atoms: &str, bonds: &str) -> Code {
    let findings = provider(&document(atoms, bonds)).expect_err("malformed rows are refused");
    findings[0].code()
}

#[test]
fn interleaved_rows_of_one_component_resolve_together() {
    let doc = document(
        "LIG C1 C N\nOTH N1 N N\nLIG O1 O N\nOTH N2 N N\n",
        "LIG C1 O1 SING N\nOTH N1 N2 DOUB N\n",
    );
    let embedded = provider(&doc).expect("well formed");
    let lig = embedded.get("LIG").expect("lookup").expect("defined");
    let names: Vec<&str> = lig.atoms.iter().map(|atom| atom.name.as_ref()).collect();
    assert_eq!(names, ["C1", "O1"]);
    assert_eq!(lig.bonds.len(), 1);
    assert_eq!(
        embedded
            .get("OTH")
            .expect("lookup")
            .expect("defined")
            .atoms
            .len(),
        2
    );
}

#[test]
fn a_component_is_lowered_once_and_shared() {
    let doc = document("LIG C1 C N\nLIG O1 O N\n", "LIG C1 O1 SING N\n");
    let embedded = provider(&doc).expect("well formed");
    let first = embedded.get("LIG").expect("lookup").expect("defined");
    let again = embedded
        .clone()
        .get("LIG")
        .expect("lookup")
        .expect("defined");
    assert!(Arc::ptr_eq(&first, &again));
}

#[test]
fn taking_the_tables_matches_copying_them() {
    let mut doc = document("LIG C1 C N\nLIG O1 O N\n", "LIG C1 O1 SING N\n");
    let copied = provider(&doc).expect("well formed");
    let taken =
        EmbeddedProvider::take_from(&mut doc, DictionaryVersion::new("t")).expect("well formed");
    assert!(
        doc.first_block()
            .and_then(|b| b.category("chem_comp_atom"))
            .is_none()
    );
    assert_eq!(
        *copied.get("LIG").expect("lookup").expect("defined"),
        *taken.get("LIG").expect("lookup").expect("defined")
    );
}

#[test]
fn malformed_rows_are_refused_at_construction_even_when_unused() {
    assert_eq!(first_code("LIG C1 Qq N\n", ""), Code::E2003);
    assert_eq!(first_code("LIG C1 C N\nLIG C1 O N\n", ""), Code::E2005);
    assert_eq!(first_code("LIG C1 C X\n", ""), Code::E2003);
    assert_eq!(
        first_code("LIG C1 C N\n", "LIG C1 Z9 SING N\n"),
        Code::E2005
    );
    assert_eq!(
        first_code("LIG C1 C N\n", "LIG C1 C1 SING N\n"),
        Code::E2005
    );
    assert_eq!(
        first_code("LIG C1 C N\nLIG O1 O N\n", "LIG C1 O1 TRIPLE N\n"),
        Code::E2003
    );
    assert_eq!(
        first_code(
            "LIG C1 C N\nLIG O1 O N\n",
            "LIG C1 O1 SING N\nLIG O1 C1 SING N\n"
        ),
        Code::E2005
    );
}

#[test]
fn a_repeated_component_id_is_refused() {
    let input = InputBuffer::from_bytes(
        b"data_s\nloop_\n_chem_comp.id\n_chem_comp.name\n_chem_comp.type\n\
LIG a NON-POLYMER\nLIG b NON-POLYMER\n\
loop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n\
_chem_comp_atom.pdbx_aromatic_flag\nLIG C1 C N\n"
            .to_vec(),
    );
    let doc = molframe_cif::parse(&input).expect("parses").0;
    let findings = provider(&doc).expect_err("duplicate definition");
    assert_eq!(findings[0].code(), Code::E2005);
}
