use super::{parse_stockholm, write_stockholm};

#[test]
fn interleaved_blocks_are_concatenated_per_name() {
    let text = "# STOCKHOLM 1.0\n\
seq1 ACGT\n\
seq2 AC-T\n\
\n\
seq1 GGCC\n\
seq2 GG-C\n\
//\n";
    let records = parse_stockholm(text);
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].id, "seq1");
    assert_eq!(records[0].sequence, b"ACGTGGCC");
    assert_eq!(records[1].sequence, b"AC-TGG-C");
}

#[test]
fn markup_lines_are_ignored() {
    let text = "# STOCKHOLM 1.0\n\
#=GC SS_cons ....\n\
seq1 ACGT\n\
//\n";
    let records = parse_stockholm(text);
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].sequence, b"ACGT");
}

#[test]
fn writing_then_reading_returns_the_same_records() {
    let original = parse_stockholm("# STOCKHOLM 1.0\nseq1 AC-GT\nseq2 ACGGT\n//\n");
    let round_tripped = parse_stockholm(&write_stockholm(&original));
    assert_eq!(round_tripped, original);
}
