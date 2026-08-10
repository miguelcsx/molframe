use super::{parse_phylip, write_phylip};

#[test]
fn the_header_is_skipped_and_each_line_is_one_record() {
    let records = parse_phylip(" 2 8\nseq1 ACGTACGT\nseq2 ACG-ACGT\n");
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].id, "seq1");
    assert_eq!(records[0].sequence, b"ACGTACGT");
    assert_eq!(records[1].sequence, b"ACG-ACGT");
}

#[test]
fn writing_then_reading_returns_the_same_records() {
    let original = parse_phylip(" 2 5\nseq1 AC-GT\nseq2 ACGGT\n");
    let round_tripped = parse_phylip(&write_phylip(&original));
    assert_eq!(round_tripped, original);
}
