use super::*;
use crate::{Provided, TableIdentity, TableKind, VdwRadii};

#[test]
fn atomic_weight_period_group_and_valence_are_defined_for_every_element() {
    let Some(carbon) = element_properties(Element::CARBON) else {
        panic!("carbon is known");
    };
    assert!((carbon.atomic_weight - 12.011).abs() < f64::EPSILON);
    for atomic_number in 1..=118 {
        let element = Element::from_atomic_number(atomic_number);
        let Some(properties) = element_properties(element) else {
            panic!("properties {atomic_number} absent")
        };
        assert!(properties.atomic_weight > 0.0);
        assert!((1..=7).contains(&properties.period));
        assert!(properties.valence_electrons > 0);
    }
}

fn radii(rows: &str) -> RadiusSource {
    let text = format!("name: lab\nversion: 1\nsource: s\nlicence: CC0-1.0\n{rows}");
    match VdwRadii::decode(text.as_bytes(), TableIdentity::supplied("x", "x", "x", "x")) {
        Ok(Provided { identity, table }) => RadiusSource::from(Provided { identity, table }),
        Err(error) => panic!("table refused: {error}"),
    }
}

#[test]
fn per_atom_radii_follow_the_table_and_mark_unknown_elements() {
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
    let Ok(record) = crate::parse_mol_record(WATER_MOL) else {
        panic!("the record parses");
    };
    let Ok(structure) = crate::mol_record_to_structure(&record) else {
        panic!("the record lowers");
    };
    let radii = atom_radii(&structure, &radii("O 1.5\nH 1.1\n"));
    assert_eq!(radii.len(), 3);
    assert!((radii[0] - 1.5).abs() < 1e-6);
    assert!((radii[1] - 1.1).abs() < 1e-6);
    assert!(radii[2].is_nan(), "an unknown element has no radius");
}
