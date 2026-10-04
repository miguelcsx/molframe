use super::*;

#[test]
fn atomic_weight_and_table_identity_are_stable() {
    let carbon = element_properties(Element::CARBON).expect("carbon is known");
    assert!((carbon.atomic_weight - 12.011).abs() < f64::EPSILON);
    assert_eq!(RadiusSet::Bondi.table().version, "1964");
}

#[test]
fn radius_sets_do_not_silently_fall_back_to_each_other() {
    assert_eq!(vdw_radius(Element::CARBON, RadiusSet::Bondi), Some(1.70));
    assert_eq!(vdw_radius(Element::IRON, RadiusSet::Bondi), None);
    assert_eq!(vdw_radius(Element::IRON, RadiusSet::Alvarez), Some(2.44));
}

#[test]
fn every_radius_set_round_trips_through_its_name() {
    for name in RadiusSet::NAMES {
        let set: RadiusSet = name.parse().expect("a listed name parses");
        assert_eq!(set.name(), *name);
    }
    assert_eq!("amber_united".parse(), Ok(RadiusSet::AmberUnited));
    let refused = "uff".parse::<RadiusSet>().expect_err("not a radius set");
    assert_eq!(refused.field, "radii");
}
