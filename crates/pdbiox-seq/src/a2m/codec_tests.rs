use super::*;

#[test]
fn insertion_columns_preserve_raw_data_but_not_match_projection() {
    let text = ">a\nACgt-D\n>b\nA-..CD\n";
    let Ok(records) = parse_a2m(text) else {
        panic!("valid A2M");
    };
    assert_eq!(records[0].sequence, b"ACgt-D");
    assert_eq!(a2m_match_columns(&records[0].sequence), b"AC-D");
    assert_eq!(a2m_match_columns(&records[1].sequence), b"A-CD");
}

#[test]
fn unequal_match_coordinates_are_rejected() {
    assert!(matches!(
        parse_a2m(">a\nACD\n>b\nAC\n"),
        Err(A2mError::MatchLength { .. })
    ));
}
