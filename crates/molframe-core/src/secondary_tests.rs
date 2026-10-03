use super::*;

#[test]
fn every_helix_class_is_a_helix_and_nothing_else_is() {
    for state in [
        SecondaryStructure::AlphaHelix,
        SecondaryStructure::ThreeTenHelix,
        SecondaryStructure::PiHelix,
        SecondaryStructure::OtherHelix,
        SecondaryStructure::PolyProline,
    ] {
        assert!(state.is_helix(), "{state:?}");
    }
    for state in [
        SecondaryStructure::Unknown,
        SecondaryStructure::Coil,
        SecondaryStructure::BetaBridge,
        SecondaryStructure::Strand,
        SecondaryStructure::Turn,
        SecondaryStructure::Bend,
    ] {
        assert!(!state.is_helix(), "{state:?}");
    }
}

#[test]
fn a_deposited_state_outranks_an_inferred_one() {
    assert!(SecondarySource::File.rank() > SecondarySource::Dssp.rank());
    assert!(SecondarySource::Dssp.rank() > SecondarySource::CaOnly.rank());
    assert!(SecondarySource::CaOnly.rank() > SecondarySource::None.rank());
}

#[test]
fn polyproline_appends_its_code_and_beta_helpers_are_distinct() {
    assert_eq!(SecondaryStructure::PolyProline.code(), 10);
    assert_eq!(SecondaryStructure::OtherHelix.code(), 7);
    assert_eq!(SecondaryStructure::BetaBridge.code(), 8);
    assert!(SecondaryStructure::Strand.is_strand());
    assert!(SecondaryStructure::Strand.is_sheet_like());
    assert!(!SecondaryStructure::BetaBridge.is_strand());
    assert!(SecondaryStructure::BetaBridge.is_sheet_like());
    assert!(!SecondaryStructure::PolyProline.is_sheet_like());
}
