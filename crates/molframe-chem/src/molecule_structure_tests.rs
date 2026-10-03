use super::*;
use crate::{parse_mol_record, parse_mol2_record};

const WATER_MOL: &str = "water\n  test\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n\
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n\
    0.7570    0.5860    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n\
   -0.7570    0.5860    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n\
  1  2  1  0\n  1  3  1  0\nM  END\n";

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
