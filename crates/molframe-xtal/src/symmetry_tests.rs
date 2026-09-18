use super::{Rational, SYMMETRY_EXTENSION, SymmetryExt, SymmetryOperation, lower_symmetry};
use crate::view::tests::{ENTRY, attached};
use molframe_cif::parse;
use molframe_core::{Code, InputBuffer};

const SYMMETRY: &str = r"data_symmetry
_space_group.IT_number 15
_space_group.name_H-M_alt 'C 2/c'
_space_group.name_Hall '-C 2yc'
_space_group.crystal_system monoclinic
loop_
_space_group_symop.id
_space_group_symop.operation_xyz
1 'x,y,z'
2 '-x,-y,-z'
3 '-x,1/2+y,1/2-z'
4 'x,1/2-y,1/2+z'
";

#[test]
fn exact_parser_handles_official_and_mixed_axis_forms() {
    let operation = match SymmetryOperation::parse("3", "-y+x,-y,1/3+z") {
        Ok(operation) => operation,
        Err(finding) => panic!("operation failed: {finding}"),
    };
    assert_eq!(operation.rotation, [[1, -1, 0], [0, -1, 0], [0, 0, 1]]);
    assert_eq!(
        operation.translation[2],
        Rational::new(1, 3).ok().unwrap_or(Rational::ZERO)
    );
    let actual = operation.apply_fractional([0.25, 0.5, 0.1]);
    assert_near(actual, [-0.25, -0.5, 0.1 + 1.0 / 3.0]);
}

#[test]
fn lowering_keeps_group_setting_and_file_order() {
    let set = parsed_set(SYMMETRY);
    assert_eq!(set.international_number, Some(15));
    assert_eq!(set.hermann_mauguin.as_deref(), Some("C 2/c"));
    assert_eq!(set.hall.as_deref(), Some("-C 2yc"));
    assert_eq!(set.operations().len(), 4);
    assert!(set.operations()[0].is_identity());
    assert_eq!(set.operations()[2].id.as_ref(), "3");
}

#[test]
fn legacy_category_and_structure_extension_are_supported() {
    let legacy = "data_old\n_symmetry.Int_Tables_number 1\n_symmetry.space_group_name_H-M 'P 1'\n_symmetry_equiv.id 1\n_symmetry_equiv.pos_as_xyz 'x,y,z'\n";
    let set = parsed_set(legacy);
    let structure = attached(ENTRY).with_extension(SYMMETRY_EXTENSION, set);
    assert_eq!(
        structure
            .symmetry_set()
            .and_then(|set| set.international_number),
        Some(1)
    );
}

#[test]
fn singular_rotations_and_invalid_group_numbers_are_diagnostics() {
    let invalid = SYMMETRY
        .replace("_space_group.IT_number 15", "_space_group.IT_number 231")
        .replace("'x,y,z'", "'x,x,z'");
    let input = InputBuffer::from_bytes(invalid.into_bytes());
    let document = match parse(&input) {
        Ok((document, _)) => document,
        Err(findings) => panic!("fixture parse failed: {findings:?}"),
    };
    let Err(findings) = lower_symmetry(&document) else {
        panic!("invalid symmetry accepted")
    };
    assert!(findings.iter().any(|finding| finding.code() == Code::E6015));
}

#[test]
fn metadata_without_explicit_operations_resolves_the_exact_hall_setting() {
    let set = parsed_set(
        "data_setting\n_space_group.IT_number 230\n_space_group.name_Hall '-I 4bd 2c 3'\n",
    );
    assert_eq!(set.hall_number, Some(530));
    assert_eq!(set.operations().len(), 96);

    let canonical = parsed_set("data_type\n_space_group.IT_number 1\n");
    assert_eq!(canonical.hall_number, Some(1));
    assert_eq!(canonical.operations().len(), 1);
    assert!(canonical.operations()[0].is_identity());
}

fn parsed_set(text: &str) -> super::SymmetrySet {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let document = match parse(&input) {
        Ok((document, _)) => document,
        Err(findings) => panic!("fixture parse failed: {findings:?}"),
    };
    match lower_symmetry(&document) {
        Ok(set) => set,
        Err(findings) => panic!("symmetry lowering failed: {findings:?}"),
    }
}

fn assert_near(actual: [f64; 3], expected: [f64; 3]) {
    assert!(
        actual
            .iter()
            .zip(expected)
            .all(|(actual, expected)| (*actual - expected).abs() < 1e-12)
    );
}
