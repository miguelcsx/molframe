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

#[test]
fn per_atom_radii_follow_the_requested_set_and_mark_unknown_elements() {
    const WATER_MOL: &str = "water
  test

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7570    0.5860    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7570    0.5860    0.0000 Xx  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
M  END
";
    let record = crate::parse_mol_record(WATER_MOL).expect("the record parses");
    let structure = crate::mol_record_to_structure(&record).expect("the record lowers");
    let radii = atom_radii(&structure, RadiusSet::Bondi);
    assert_eq!(radii.len(), 3);
    assert_eq!(
        Some(radii[0]),
        vdw_radius(Element::OXYGEN, RadiusSet::Bondi)
    );
    assert_eq!(
        Some(radii[1]),
        vdw_radius(Element::HYDROGEN, RadiusSet::Bondi)
    );
    assert!(radii[2].is_nan(), "an unknown element has no radius");
}
