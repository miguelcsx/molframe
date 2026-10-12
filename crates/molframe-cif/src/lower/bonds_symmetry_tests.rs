use super::{OperatorTable, is_local, is_stated, parse, unrepresentable};
use crate::document::DataBlock;

fn block(text: &str) -> DataBlock {
    let input = molframe_core::io::InputBuffer::from_bytes(text.as_bytes().to_vec());
    let (document, findings) = crate::parse(&input).unwrap_or_else(|error| panic!("{error:?}"));
    assert!(findings.is_empty(), "{findings:?}");
    document
        .first_block()
        .cloned()
        .unwrap_or_else(|| panic!("one block"))
}

const IDENTITY_TABLE: &str = "data_t\nloop_\n_space_group_symop.id\n\
_space_group_symop.operation_xyz\n1 'x,y,z'\n";

#[test]
fn an_offset_digit_is_a_translation_not_an_operator_digit() {
    // `_555` is no translation; the middle digit counts along b, the last along c.
    assert_eq!(parse("1_555"), Some((1, [0, 0, 0])));
    assert_eq!(parse("1_565"), Some((1, [0, 1, 0])));
    assert_eq!(parse("1_556"), Some((1, [0, 0, 1])));
    assert_eq!(parse("2_555"), Some((2, [0, 0, 0])));
    assert_eq!(parse("1_444"), Some((1, [-1, -1, -1])));
    assert_eq!(parse("1"), Some((1, [0, 0, 0])));
    assert_eq!(parse("n_555"), None);
    assert_eq!(parse("1_55"), None);
    assert_eq!(parse("1_5a5"), None);
}

#[test]
fn unknown_and_inapplicable_are_not_stated_values() {
    assert!(!is_stated("?"));
    assert!(!is_stated("."));
    assert!(!is_stated(""));
    assert!(!is_stated("   "));
    assert!(is_stated("1_555"));
}

#[test]
fn the_local_identity_is_recognized_with_and_without_a_table() {
    let empty = OperatorTable::of(&block("data_t\n_entry.id x\n"));
    assert!(is_local("1", &empty));
    assert!(is_local("1_555", &empty));
    assert!(!is_local("2_555", &empty));
    assert!(!is_local("1_565", &empty));
    assert!(!is_local("garbage", &empty));

    let identity = OperatorTable::of(&block(IDENTITY_TABLE));
    assert!(is_local("1", &identity));
    assert!(is_local("1_555", &identity));
    assert!(!is_local("2", &identity));
}

#[test]
fn a_table_that_defines_operator_one_as_a_translation_is_honoured() {
    // Operator 1 is not identity here, so a connection naming it is not local
    // and cannot be lowered to a bond between the atoms of this asymmetric unit.
    let shifted = block(
        "data_t\nloop_\n_space_group_symop.id\n_space_group_symop.operation_xyz\n\
1 'x,y+1/2,z'\n2 '-x,-y,-z'\n",
    );
    let table = OperatorTable::of(&shifted);
    assert!(!is_local("1", &table));
    assert!(!is_local("1_555", &table));
    assert!(!is_local("2", &table));
}

#[test]
fn whitespace_and_leading_plus_do_not_hide_the_identity() {
    let table = OperatorTable::of(&block(
        "data_t\nloop_\n_space_group_symop.id\n_space_group_symop.operation_xyz\n\
1 ' x, y , z '\n",
    ));
    assert!(is_local("1_555", &table));
    let legacy = OperatorTable::of(&block(
        "data_t\nloop_\n_symmetry_equiv.id\n_symmetry_equiv.pos_as_xyz\n1 'X,Y,Z'\n",
    ));
    assert!(is_local("1_555", &legacy));
}
#[test]
fn a_nonzero_translation_on_either_partner_is_refused() {
    let table = OperatorTable::of(&block(IDENTITY_TABLE));
    let found = unrepresentable(Some("1_555"), Some("1_565"), &table)
        .unwrap_or_else(|| panic!("the second partner is on the next cell along b"));
    assert_eq!(found.partner, "ptnr2_symmetry");
    assert_eq!(found.value, "1_565");
    let found = unrepresentable(Some("2_555"), Some("1_555"), &table)
        .unwrap_or_else(|| panic!("operator two is not the identity"));
    assert_eq!(found.partner, "ptnr1_symmetry");
    assert!(unrepresentable(Some("1_555"), Some("1_555"), &table).is_none());
    assert!(unrepresentable(None, Some("?"), &table).is_none());
    assert!(unrepresentable(Some("."), None, &table).is_none());
    let malformed = unrepresentable(Some("1_55"), None, &table)
        .unwrap_or_else(|| panic!("a code that cannot be parsed is not a local claim"));
    assert_eq!(malformed.value, "1_55");
}
