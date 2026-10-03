use super::*;

#[test]
fn every_helix_class_is_a_helix_and_nothing_else_is() {
    for state in [
        SecondaryStructure::AlphaHelix,
        SecondaryStructure::ThreeTenHelix,
        SecondaryStructure::PiHelix,
        SecondaryStructure::OtherHelix,
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
