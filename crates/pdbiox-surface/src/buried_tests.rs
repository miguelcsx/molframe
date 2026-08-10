use super::{MoleculeRole, buried_solvent_excluded_surface, buried_surface};

#[test]
fn groups_pulled_far_apart_bury_nothing() {
    let positions = [[0.0, 0.0, 0.0], [100.0, 0.0, 0.0]];
    let radii = [1.5, 1.5];
    let Ok(result) = buried_surface(&positions, &radii, 1.4, 400, &[true, false]) else {
        panic!("valid");
    };
    assert!(result.buried.abs() < 1e-6, "buried {}", result.buried);
    assert!((result.together - result.first_alone - result.second_alone).abs() < 1e-6);
}

#[test]
fn overlapping_groups_bury_a_positive_area() {
    // Two atoms whose expanded spheres overlap: forming the pair hides area.
    let positions = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let radii = [1.5, 1.5];
    let Ok(result) = buried_surface(&positions, &radii, 0.5, 2000, &[true, false]) else {
        panic!("valid");
    };
    assert!(result.buried > 0.0, "expected contact to bury area");
    assert!(result.together < result.first_alone + result.second_alone);
}

#[test]
fn an_empty_partner_buries_no_area() {
    let positions = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let radii = [1.5, 1.5];
    // Everything is in the first group; the second is empty.
    let Ok(result) = buried_surface(&positions, &radii, 1.4, 400, &[true, true]) else {
        panic!("valid");
    };
    assert!(result.second_alone.abs() < 1e-12);
    assert!(result.buried.abs() < 1e-6, "buried {}", result.buried);
}

#[test]
fn a_mismatched_partition_is_rejected() {
    let positions = [[0.0, 0.0, 0.0]];
    let radii = [1.5];
    assert!(buried_surface(&positions, &radii, 1.4, 100, &[true, false]).is_err());
}

#[test]
fn ses_roles_exclude_environment_atoms_explicitly() {
    let positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [1.5, 0.0, 0.0]];
    let radii = [1.5, 1.5, 8.0];
    let roles = [
        MoleculeRole::First,
        MoleculeRole::Second,
        MoleculeRole::Excluded,
    ];
    let Ok(with_excluded) = buried_solvent_excluded_surface(&positions, &radii, &roles, 0.5, 0.3)
    else {
        panic!("explicit SES partition should be valid");
    };
    let Ok(reference) =
        buried_solvent_excluded_surface(&positions[..2], &radii[..2], &roles[..2], 0.5, 0.3)
    else {
        panic!("reference SES partition should be valid");
    };

    assert!((with_excluded.buried - reference.buried).abs() < 1e-8);
}

#[test]
fn ses_roles_are_atom_aligned() {
    let result = buried_solvent_excluded_surface(&[[0.0, 0.0, 0.0]], &[1.5], &[], 0.5, 0.3);
    assert!(result.is_err());
}
