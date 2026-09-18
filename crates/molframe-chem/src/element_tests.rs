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
