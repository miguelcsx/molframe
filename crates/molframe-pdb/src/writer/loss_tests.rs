use super::*;
use crate::read;
use molframe_core::io::{InputBuffer, ReadOptions};

fn parse(text: &str) -> Structure {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    match read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("read failed: {findings:?}"),
    }
}

#[test]
fn a_file_the_format_holds_exactly_reports_no_loss() {
    let source = "\
ATOM      1  N   GLY A   1      27.340  24.430   2.614  1.00 10.00           N
END
";
    assert!(projection_loss(&parse(source), &PdbOptions::new()).is_empty());
}

#[test]
fn a_value_is_rounded_exactly_when_the_columns_cannot_hold_it() {
    assert!(rounded(1.2345, 1e3));
    assert!(!rounded(1.234, 1e3));
    assert!(rounded(0.125, 1e2));
    assert!(!rounded(0.12, 1e2));
}
