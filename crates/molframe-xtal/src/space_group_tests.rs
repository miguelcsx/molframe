use super::{
    settings, space_group_by_hall, space_group_by_hermann_mauguin_in_cell, space_group_setting,
    space_group_settings,
};
use crate::SymmetrySet;
use molframe_core::structure::UnitCell;
use std::collections::BTreeSet;

#[test]
fn catalogue_covers_all_types_and_all_530_hall_settings() {
    let catalogue = match settings() {
        Ok(catalogue) => catalogue,
        Err(finding) => panic!("catalogue failed: {finding}"),
    };
    assert_eq!(catalogue.len(), 530);
    assert_eq!(
        catalogue
            .iter()
            .map(|setting| setting.international_number)
            .collect::<BTreeSet<_>>()
            .len(),
        230
    );
    assert!(
        catalogue
            .iter()
            .all(|setting| !setting.operations.is_empty())
    );
}

#[test]
fn hall_number_symbol_and_type_lookups_preserve_setting_identity() {
    let last = match space_group_setting(530) {
        Ok(setting) => setting,
        Err(finding) => panic!("lookup failed: {finding}"),
    };
    assert_eq!(last.international_number, 230);
    assert_eq!(last.hall_symbol.as_ref(), "-I 4bd 2c 3");
    assert_eq!(last.operations.len(), 96);
    assert_eq!(
        space_group_by_hall("  -I   4bd 2c 3 ").map(|setting| setting.hall_number),
        Ok(530)
    );
    assert!(space_group_settings(3).is_ok_and(|settings| settings.len() > 1));
}

#[test]
fn every_operation_has_an_exact_inverse_representative() {
    let catalogue = match settings() {
        Ok(catalogue) => catalogue,
        Err(finding) => panic!("catalogue failed: {finding}"),
    };
    for setting in catalogue {
        let set = setting.symmetry_set();
        for operation in 0..set.operations().len() {
            assert!(
                set.inverse_image(operation, [0; 3]).is_ok(),
                "Hall setting {} operation {} has no inverse",
                setting.hall_number,
                operation
            );
        }
    }
}

#[test]
fn hermann_mauguin_lookup_accepts_pdb_and_catalogue_spellings() {
    for (spelling, hall_number_of_type) in [
        ("P 21 21 21", 19),
        ("p212121", 19),
        ("P 2_1 2_1 2_1", 19),
        ("I 2 2 2", 23),
        ("P 21 3", 198),
        ("P 1", 1),
    ] {
        let Ok(setting) = super::space_group_by_hermann_mauguin(spelling) else {
            panic!("{spelling} resolves")
        };
        assert_eq!(
            setting.international_number, hall_number_of_type,
            "{spelling}"
        );
    }
    // The full symbol names the unique axis; the short one means the first setting.
    let Ok(unique_c) = super::space_group_by_hermann_mauguin("P 1 1 21") else {
        panic!("the c-axis setting resolves")
    };
    assert_eq!(&*unique_c.hall_symbol, "P 2c");
    assert!(super::space_group_by_hermann_mauguin("Q 9 9 9").is_err());
}

#[test]
fn hall_change_of_basis_resolves_to_the_catalogue_setting() {
    for (symbol, hall_number) in [
        ("P 61 2 (x,y,z+5/12)", 472),
        ("P 61 2 (0 0 5)", 472),
        ("P 65 2 (x,y,z+1/12)", 473),
        ("P 2 2ab (y,z,x)", 114),
    ] {
        assert_eq!(
            space_group_by_hall(symbol).map(|setting| setting.hall_number),
            Ok(hall_number),
            "{symbol}"
        );
    }
    assert!(space_group_by_hall("P 2 2ab (y,z").is_err());
    assert!(space_group_by_hall("Q 9 (0 0 1)").is_err());
}

#[test]
fn hermann_mauguin_choices_and_pdb_rhombohedral_spellings() {
    let hall = |symbol: &str| {
        super::space_group_by_hermann_mauguin(symbol).map(|setting| setting.hall_number)
    };
    assert_eq!(hall("I 1 2 1"), Ok(11));
    assert_eq!(hall("H 3"), Ok(433));
    assert_eq!(hall("H 3 2"), Ok(444));
    assert_eq!(hall("R 3 :R"), Ok(434));
    assert_eq!(hall("R 3 :H"), Ok(433));
    let rhombohedral = UnitCell {
        lengths: [50.0; 3],
        angles: [80.0; 3],
    };
    let hexagonal = UnitCell {
        lengths: [50.0, 50.0, 120.0],
        angles: [90.0, 90.0, 120.0],
    };
    let in_cell = |symbol: &str, cell: &UnitCell| {
        space_group_by_hermann_mauguin_in_cell(symbol, Some(cell))
            .map(|setting| setting.hall_number)
    };
    assert_eq!(in_cell("R 3 2", &rhombohedral), Ok(445));
    assert_eq!(in_cell("R 3 2", &hexagonal), Ok(444));
    assert_eq!(in_cell("H 3 2", &rhombohedral), Ok(444));
}

#[test]
fn metadata_resolution_uses_each_symbol_for_what_it_fixes() {
    let resolve = |number: Option<u16>, hermann_mauguin: Option<&str>, hall: Option<&str>| {
        let mut set = SymmetrySet {
            international_number: number,
            hermann_mauguin: hermann_mauguin.map(Into::into),
            hall: hall.map(Into::into),
            ..SymmetrySet::default()
        };
        super::resolve_metadata_operations(&mut set, None).map(|()| set)
    };
    // The number alone would mean C 1 2 1; the symbol says the body-centred cell.
    let body_centred = resolve(Some(5), Some("I 1 2 1"), None);
    assert!(body_centred.is_ok_and(|set| set.hall_number == Some(11)));
    let screw = resolve(Some(178), None, Some(" P 61 2 (x,y,z+5/12)"));
    assert!(screw.is_ok_and(|set| set.hall_number == Some(472) && set.operations().len() == 12));
    // A number that contradicts the symbol is reported, not overridden.
    assert!(resolve(Some(4), Some("I 1 2 1"), None).is_err());
    // An unknown spelling of a type with several settings stays unresolved.
    assert!(resolve(Some(5), Some("C 2 weird"), None).is_ok_and(|set| set.operations().is_empty()));
    assert!(resolve(Some(19), Some("unknown"), None).is_ok_and(|set| set.hall_number == Some(115)));
}
