use crate::PdbHeaders;
use crate::PdbHeadersExt;
use molframe_bench::{Sample, input};
use molframe_core::io::ReadOptions;

fn headers_of(text: &str) -> PdbHeaders {
    let mut headers = PdbHeaders::default();
    for line in text.lines() {
        headers.push(line.get(..6).unwrap_or("").trim(), line);
    }
    headers
}

#[test]
fn hemoglobin_carries_one_given_operator_near_a_proper_rotation() {
    let bytes = Sample::Medium.pdb().expect("the sample ships a PDB");
    let (structure, _) = crate::read(&input(bytes), &ReadOptions::new()).expect("sample reads");
    let operators = structure.pdb_headers().expect("headers kept").mtrix();
    assert_eq!(operators.len(), 1);
    let operator = &operators[0];
    assert!(operator.given);
    assert!((operator.rotation[0][0] + 0.952_899).abs() < 1e-9);
    assert!((operator.translation[2] - 81.782_29).abs() < 1e-9);
}

#[test]
fn an_operator_missing_a_row_is_dropped_and_a_blank_flag_is_generated() {
    let text = "\
MTRIX1   1  1.000000  0.000000  0.000000        1.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        2.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        3.00000    1
MTRIX1   2  1.000000  0.000000  0.000000        1.00000
MTRIX2   2  0.000000  1.000000  0.000000        2.00000
MTRIX3   2  0.000000  0.000000  1.000000        3.00000
MTRIX1   3  1.000000  0.000000  0.000000        1.00000
MTRIX2   3  0.000000  1.000000  0.000000        2.00000
";
    let operators = headers_of(text).mtrix();
    assert_eq!(operators.len(), 2);
    assert!(operators[0].given);
    assert!(!operators[1].given);
    assert!(
        operators[0]
            .translation
            .iter()
            .zip([1.0, 2.0, 3.0])
            .all(|(actual, expected)| (actual - expected).abs() < 1e-12)
    );
}

#[test]
fn garbage_columns_yield_no_operator() {
    assert!(
        headers_of("MTRIX1   x  abc\nMTRIX2\nMTRIX3   1")
            .mtrix()
            .is_empty()
    );
}
