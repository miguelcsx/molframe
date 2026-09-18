use super::{a3m_match_columns, parse_a3m};

#[test]
fn insertions_and_their_dots_are_dropped_from_the_match_columns() {
    // Upper-case A,C,-,D,E,F are match columns; lower-case g,t and the dot are
    // insertions relative to the query.
    assert_eq!(a3m_match_columns(b"AC-gtDE.F"), b"AC-DEF");
}

#[test]
fn two_records_reduce_to_the_same_match_length() {
    let text = ">query\nACDEF\n>hit\nACdeDEF\n";
    let Ok(records) = parse_a3m(text) else {
        panic!("valid A3M");
    };
    assert_eq!(records.len(), 2);
    let query = a3m_match_columns(&records[0].sequence);
    let hit = a3m_match_columns(&records[1].sequence);
    assert_eq!(query.len(), hit.len());
    assert_eq!(query, b"ACDEF");
    assert_eq!(hit, b"ACDEF");
}
