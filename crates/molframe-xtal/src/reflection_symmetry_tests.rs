use super::*;
use crate::{Rational, space_group_by_hall};

#[test]
fn screw_glide_and_centering_extinctions_preserve_allowed_reflections() {
    let screw = space_group_by_hall("P 2yb").expect("P21").symmetry_set();
    assert!(
        screw
            .reflection_symmetry([0, 1, 0])
            .expect("valid")
            .systematically_absent
    );
    assert!(
        !screw
            .reflection_symmetry([0, 2, 0])
            .expect("valid")
            .systematically_absent
    );
    assert!(
        !screw
            .reflection_symmetry([1, 1, 0])
            .expect("valid")
            .systematically_absent
    );
    let centered = space_group_by_hall("I 2 2 3").expect("I23").symmetry_set();
    assert!(
        centered
            .reflection_symmetry([1, 0, 0])
            .expect("valid")
            .systematically_absent
    );
    assert!(
        !centered
            .reflection_symmetry([1, 1, 0])
            .expect("valid")
            .systematically_absent
    );
    let glide = space_group_by_hall("P -2yc").expect("Pc").symmetry_set();
    assert!(
        glide
            .reflection_symmetry([1, 0, 1])
            .expect("valid")
            .systematically_absent
    );
    assert!(
        !glide
            .reflection_symmetry([1, 1, 1])
            .expect("valid")
            .systematically_absent
    );
}

#[test]
fn centricity_and_epsilon_include_complete_group_centering() {
    let inversion = space_group_by_hall("-P 1")
        .expect("inversion")
        .symmetry_set();
    assert_eq!(
        inversion.reflection_symmetry([1, 2, 3]).expect("valid"),
        ReflectionSymmetry {
            centric: true,
            systematically_absent: false,
            epsilon_factor: 1,
        }
    );
    let centered = space_group_by_hall("I 2 2 3").expect("I23").symmetry_set();
    let origin = centered.reflection_symmetry([0, 0, 0]).expect("origin");
    assert!(origin.centric);
    assert!(!origin.systematically_absent);
    assert_eq!(origin.epsilon_factor, centered.operations().len());
    assert_eq!(
        centered
            .reflection_symmetry([1, 2, 3])
            .expect("valid")
            .epsilon_factor,
        2
    );
    assert!(
        SymmetrySet::default()
            .reflection_symmetry([1, 0, 0])
            .is_err()
    );
}

#[test]
fn reciprocal_rotations_use_the_transpose_and_report_overflow() {
    let operation = SymmetryOperation::parse("1", "y,z,x").expect("rotation");
    assert_eq!(operation.apply_to_hkl([2, 3, 5]).expect("valid"), [5, 2, 3]);
    let inversion = SymmetryOperation::parse("2", "-x,-y,-z").expect("inversion");
    assert!(inversion.apply_to_hkl([i32::MIN, 0, 0]).is_err());
}

#[test]
fn translation_phases_are_exact_for_large_signed_indices_and_denominators() {
    let operation = SymmetryOperation {
        id: "1".into(),
        rotation: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        translation: [
            Rational::new(1, 2).expect("half"),
            Rational::new(1, 3).expect("third"),
            Rational::new(1, 6).expect("sixth"),
        ],
    };
    assert!(!operation.has_reflection_phase_shift([i32::MAX, i32::MAX, i32::MAX]));
    assert!(!operation.has_reflection_phase_shift([i32::MIN, i32::MIN, i32::MIN]));
    assert!(operation.has_reflection_phase_shift([-1, 0, 0]));
    let operation = SymmetryOperation {
        translation: [
            Rational::new(i32::MAX, u32::MAX).expect("fraction"),
            Rational::new(i32::MAX, u32::MAX - 1).expect("fraction"),
            Rational::new(i32::MAX, u32::MAX - 2).expect("fraction"),
        ],
        ..operation
    };
    assert!(operation.has_reflection_phase_shift([i32::MIN, 0, 0]));
}
