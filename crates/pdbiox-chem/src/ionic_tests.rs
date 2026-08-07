use super::{IonicSpin, ionic_radii};
use crate::element_properties;
use pdbiox_core::Element;

#[test]
fn all_118_elements_have_complete_identity_and_periodic_properties() {
    for atomic_number in 1..=118 {
        let element = Element::from_atomic_number(atomic_number);
        let Some(properties) = element_properties(element) else {
            panic!("properties {atomic_number} absent")
        };
        assert!(properties.atomic_weight > 0.0);
        assert!(
            properties
                .covalent_radius
                .is_some_and(|radius| radius > 0.0)
        );
        assert!((1..=7).contains(&properties.period));
        assert!(properties.valence_electrons > 0);
    }
}

#[test]
fn ionic_entries_preserve_charge_coordination_spin_and_absence() {
    let iron = match ionic_radii(Element::IRON) {
        Ok(entries) => entries,
        Err(finding) => panic!("ionic table failed: {finding}"),
    };
    assert!(iron.iter().any(|entry| {
        entry.charge == 2
            && entry.coordination.as_ref() == "VI"
            && entry.spin == IonicSpin::High
            && entry.ionic_radius > 0.0
    }));
    let oganesson = ionic_radii(Element::from_atomic_number(118));
    assert!(oganesson.is_ok_and(<[super::IonicRadius]>::is_empty));
}
