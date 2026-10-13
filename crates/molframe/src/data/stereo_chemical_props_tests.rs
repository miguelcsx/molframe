use super::*;

const SAMPLE: &str = "Bond\t\t\tResidue\t\tMean\t\tStdDev\n\
CA-CB\t\t\tALA\t\t1.520\t\t0.021\n\
-\n\n\
Angle\t\t\tResidue\t\tMean\t\tStdDev\n\
N-CA-CB\t\t\tALA\t\t110.1\t\t1.4\n\
-\n\n\
Non-bonded distance     Minimum Dist    Tolerance\n\
C-C                     3.4             1.5\n\
-\n";

#[test]
fn sections_become_typed_rows_and_the_non_bonded_block_is_skipped() {
    let rows = rows(SAMPLE).expect("rows");
    assert_eq!(
        rows,
        [
            "bond\tALA\tCA-CB\t1.520\t0.021",
            "angle\tALA\tN-CA-CB\t110.1\t1.4"
        ]
    );
}

#[test]
fn a_malformed_row_names_its_line() {
    let bad = SAMPLE.replace("1.520", "abc");
    let message = rows(&bad).expect_err("refused");
    assert!(
        message.contains("line 2") && message.contains("abc"),
        "{message}"
    );
    let wrong_arity = SAMPLE.replace("N-CA-CB", "N-CA");
    assert!(rows(&wrong_arity).expect_err("arity").contains("line 6"));
    assert!(rows("nothing here\n").is_err());
}
