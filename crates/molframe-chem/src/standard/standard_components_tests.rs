use super::standard_component_kind;
use crate::ComponentKind;

#[test]
fn standard_polymer_residues_are_recognised_by_name() {
    assert_eq!(standard_component_kind("ALA"), ComponentKind::AminoAcid);
    assert_eq!(standard_component_kind("DA"), ComponentKind::Nucleotide);
}

#[test]
fn curated_ccd_sugars_share_the_carbohydrate_inventory_classification() {
    for name in ["NAG", "MAN", "BMA", "FUC", "GAL"] {
        assert_eq!(standard_component_kind(name), ComponentKind::Saccharide);
    }
    assert_eq!(
        standard_component_kind("UNRECOGNIZED"),
        ComponentKind::Unknown
    );
}

#[test]
fn water_ions_and_ligands_are_left_to_the_files_declared_entities() {
    for name in ["HOH", "ZN", "ATP", "HEM"] {
        assert_eq!(
            standard_component_kind(name),
            ComponentKind::Unknown,
            "{name}"
        );
    }
}
