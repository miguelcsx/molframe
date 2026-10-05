use super::*;
use crate::{parse_mol_record, parse_mol2_record};

const WATER_MOL: &str = "water
  test

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7570    0.5860    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7570    0.5860    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
M  END
";

const WATER_MOL2: &str = "@<TRIPOS>MOLECULE\nwat\n 3 2 1 0 0\nSMALL\nNO_CHARGES\n\n\n\
@<TRIPOS>ATOM\n\
 1 OW 0.0 0.0 0.0 O.3 1 WAT 0.0\n\
 2 HW1 0.757 0.586 0.0 H 1 WAT 0.0\n\
 3 HW2 -0.757 0.586 0.0 H 1 WAT 0.0\n\
@<TRIPOS>BOND\n 1 1 2 1\n 2 1 3 ar\n";

#[test]
fn a_mol_record_becomes_one_residue_with_its_bonds() {
    let record = parse_mol_record(WATER_MOL).expect("water parses");
    let structure = mol_record_to_structure(&record).expect("water lowers");
    assert_eq!(structure.atom_count(), 3);
    assert_eq!(structure.data().bonds.len(), 2);
    let residues: Vec<_> = structure.data().residues().collect();
    assert_eq!(residues.len(), 1);
    assert_eq!(residues[0].name(), Some("WAT"));
    let names: Vec<_> = structure
        .data()
        .atoms()
        .filter_map(molframe_core::structure::AtomRef::name)
        .collect();
    assert_eq!(names, ["O1", "H1", "H2"]);
}

#[test]
fn a_mol2_record_keeps_atom_names_and_bond_types() {
    let record = parse_mol2_record(WATER_MOL2).expect("water parses");
    let structure = mol2_record_to_structure(&record).expect("water lowers");
    let names: Vec<_> = structure
        .data()
        .atoms()
        .filter_map(molframe_core::structure::AtomRef::name)
        .collect();
    assert_eq!(names, ["OW", "HW1", "HW2"]);
    let orders: Vec<_> = structure
        .data()
        .bonds
        .iter()
        .map(|bond| bond.order)
        .collect();
    assert_eq!(orders, [BondOrder::Single, BondOrder::Aromatic]);
}

#[test]
fn a_molecule_round_trips_through_a_structure() {
    let record = parse_mol_record(WATER_MOL).expect("water parses");
    let structure = molecule_to_structure("water", &record.molecule).expect("lowers");
    assert_eq!(
        structure_to_molecule(&structure).expect("extracts"),
        record.molecule
    );
}

#[test]
fn a_bond_to_an_absent_atom_is_refused() {
    let molecule = Molecule {
        atoms: vec![MolAtom {
            element: Element::UNKNOWN,
            position: [0.0; 3],
        }],
        bonds: vec![MolBond {
            first: 0,
            second: 3,
            order: 1,
        }],
    };
    assert!(molecule_to_structure("x", &molecule).is_err());
}

#[test]
fn an_unknown_bond_order_is_not_written_as_single() {
    let mut molecule = parse_mol_record(WATER_MOL).expect("parses").molecule;
    molecule.bonds[0].order = 9;
    let structure = molecule_to_structure("w", &molecule).expect("lowers");
    assert!(structure_to_molecule(&structure).is_err());
}

#[test]
fn sdf_formal_charges_and_annotations_agree_including_neutral_atoms() {
    use molframe_core::{AtomAnnotation, FORMAL_CHARGE_ANNOTATION};
    let mut record = parse_mol_record(WATER_MOL).expect("parses");
    record.atom_metadata[1].formal_charge = Some(1);
    record.atom_metadata[2].formal_charge = Some(-1);
    let structure = mol_record_to_structure(&record).expect("lowers");
    let AtomAnnotation::Integer(column) = structure
        .data()
        .annotations
        .get(FORMAL_CHARGE_ANNOTATION)
        .expect("charge annotation")
    else {
        panic!("integer charge column");
    };
    assert_eq!(
        structure
            .data()
            .atoms()
            .map(molframe_core::structure::AtomRef::formal_charge)
            .collect::<Vec<_>>(),
        vec![Some(0), Some(1), Some(-1)]
    );
    for (index, charge) in [0, 1, -1].into_iter().enumerate() {
        assert_eq!(
            column.get(u32::try_from(index).expect("small index")),
            Some((charge, Presence::Present))
        );
    }
    for (pattern, expected) in [("[-1]", 2), ("[+1]", 1), ("[O;+0]", 0)] {
        let matches = crate::SmartsPattern::parse(pattern)
            .expect("pattern")
            .find_structure_matches(&structure)
            .expect("matches");
        assert_eq!(matches[0].atom_indices.as_ref(), &[expected]);
        assert_eq!(matches.len(), 1);
    }
}

#[test]
fn benzene_fixtures_match_aromatic_atoms_and_bonds_without_rewriting_bonds() {
    for text in [
        include_str!("../tests/fixtures/benzene-aromatic.sdf"),
        include_str!("../tests/fixtures/benzene-kekulized.sdf"),
    ] {
        let record = crate::parse_sdf_records(text).expect("parses").remove(0);
        let structure = mol_record_to_structure(&record).expect("lowers");
        for pattern in ["[a]", "[c]"] {
            assert_eq!(
                crate::SmartsPattern::parse(pattern)
                    .expect("pattern")
                    .find_structure_matches(&structure)
                    .expect("matches")
                    .len(),
                6
            );
        }
        assert_eq!(
            crate::SmartsPattern::parse("c:c")
                .expect("pattern")
                .find_structure_matches(&structure)
                .expect("matches")
                .len(),
            12
        );
        assert_eq!(
            structure_to_molecule(&structure).expect("graph"),
            record.molecule
        );
    }
}

#[test]
fn a_kekulized_pyrrole_matches_its_implicit_nitrogen_hydrogen() {
    let record = crate::parse_sdf_records(include_str!("../tests/fixtures/pyrrole-kekulized.sdf"))
        .expect("parses")
        .remove(0);
    let structure = mol_record_to_structure(&record).expect("lowers");
    let matches = crate::SmartsPattern::parse("[nH]")
        .expect("pattern")
        .find_structure_matches(&structure)
        .expect("matches");
    assert_eq!(matches.len(), 1);
    assert_eq!(matches[0].atom_indices.as_ref(), &[0]);
}

#[test]
fn unproved_aromaticity_remains_unknown_and_missing_metadata_is_refused() {
    use molframe_core::{AROMATIC_ATOM_ANNOTATION, AtomAnnotation};
    let mut record = parse_mol_record(WATER_MOL).expect("parses");
    let structure = mol_record_to_structure(&record).expect("lowers");
    let AtomAnnotation::Boolean(column) = structure
        .data()
        .annotations
        .get(AROMATIC_ATOM_ANNOTATION)
        .expect("aromaticity annotation")
    else {
        panic!("boolean aromaticity column");
    };
    assert_eq!(column.get(0), Some((false, Presence::Unknown)));
    record.atom_metadata.pop();
    assert!(mol_record_to_structure(&record).is_err());
}

#[test]
fn v3000_omitted_charge_is_known_neutral_in_the_structure() {
    let text = "neutral\nmolframe\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n\
M  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\n\
M  V30 1 O 0 0 0 0\nM  V30 END ATOM\nM  V30 END CTAB\nM  END\n";
    let record = parse_mol_record(text).expect("parses");
    let structure = mol_record_to_structure(&record).expect("lowers");
    assert_eq!(
        structure
            .data()
            .atoms()
            .next()
            .expect("atom")
            .formal_charge(),
        Some(0)
    );
    assert_eq!(
        crate::SmartsPattern::parse("[O;+0]")
            .expect("pattern")
            .find_structure_matches(&structure)
            .expect("matches")
            .len(),
        1
    );
}
