use super::*;

#[test]
fn a_backbone_carbonyl_is_double_in_either_argument_order() {
    assert_eq!(standard_bond_order("ALA", "C", "O"), BondOrder::Double);
    assert_eq!(standard_bond_order("ALA", "O", "C"), BondOrder::Double);
}

#[test]
fn a_backbone_carbonyl_outside_an_amino_acid_is_single() {
    assert_eq!(standard_bond_order("LIG", "C", "O"), BondOrder::Single);
}

#[test]
fn ring_bonds_of_aromatic_side_chains_are_aromatic() {
    assert_eq!(standard_bond_order("PHE", "CG", "CD1"), BondOrder::Aromatic);
    assert_eq!(
        standard_bond_order("TRP", "CD2", "CE2"),
        BondOrder::Aromatic
    );
    assert_eq!(
        standard_bond_order("HIS", "ND1", "CE1"),
        BondOrder::Aromatic
    );
}

#[test]
fn a_ring_to_substituent_bond_is_not_aromatic() {
    assert_eq!(standard_bond_order("TYR", "CZ", "OH"), BondOrder::Single);
    assert_eq!(standard_bond_order("PHE", "CB", "CG"), BondOrder::Single);
}

#[test]
fn nucleobase_rings_are_aromatic_and_their_carbonyls_double() {
    assert_eq!(standard_bond_order("DG", "N9", "C8"), BondOrder::Aromatic);
    assert_eq!(standard_bond_order("U", "C4", "O4"), BondOrder::Double);
    assert_eq!(standard_bond_order("DA", "OP1", "P"), BondOrder::Double);
}

#[test]
fn a_bond_the_table_does_not_name_is_single() {
    assert_eq!(standard_bond_order("LYS", "CE", "NZ"), BondOrder::Single);
}

#[test]
fn common_water_names_are_recognised() {
    assert!(is_water_component("HOH"));
    assert!(is_water_component("WAT"));
    assert!(!is_water_component("HEM"));
}
