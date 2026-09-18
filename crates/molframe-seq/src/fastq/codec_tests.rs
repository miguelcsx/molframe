use super::*;

#[test]
fn wrapped_fastq_round_trips_canonically() {
    let input = "@read description\nAC\nGT\n+read\n!!\n##\n";
    let Ok(records) = parse_fastq(input) else {
        panic!("valid wrapped FASTQ");
    };
    assert_eq!(records[0].sequence, b"ACGT");
    assert_eq!(records[0].quality, b"!!##");
    let Ok(written) = write_fastq(&records) else {
        panic!("valid FASTQ records");
    };
    assert_eq!(parse_fastq(&written), Ok(records));
}

#[test]
fn truncated_quality_is_an_explicit_error() {
    assert!(matches!(
        parse_fastq("@r\nACGT\n+\n!!!\n"),
        Err(FastqError::LengthMismatch { .. })
    ));
}

#[test]
fn separator_identifier_must_match() {
    assert!(matches!(
        parse_fastq("@first\nA\n+second\n!\n"),
        Err(FastqError::MissingSeparator { .. })
    ));
}
