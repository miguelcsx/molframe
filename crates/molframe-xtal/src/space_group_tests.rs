use super::{settings, space_group_by_hall, space_group_setting, space_group_settings};
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
