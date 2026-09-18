use super::{parse_clustal, write_clustal};

#[test]
fn blocks_are_concatenated_and_conservation_lines_ignored() {
    // Built with concat! so the conservation lines keep their leading spaces —
    // a `\`-continuation would eat them and turn `** *` into a data line.
    let text = concat!(
        "CLUSTAL W (1.83) multiple sequence alignment\n",
        "\n",
        "seq1 ACGT 4\n",
        "seq2 AC-T 3\n",
        "     ** *\n",
        "\n",
        "seq1 GGCC 8\n",
        "seq2 GG-C 6\n",
        "     ** *\n",
    );
    let records = parse_clustal(text);
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].id, "seq1");
    assert_eq!(records[0].sequence, b"ACGTGGCC");
    assert_eq!(records[1].sequence, b"AC-TGG-C");
}

#[test]
fn writing_then_reading_returns_the_same_records() {
    let original = parse_clustal("CLUSTAL W\n\nseq1 AC-GT\nseq2 ACGGT\n");
    let round_tripped = parse_clustal(&write_clustal(&original));
    assert_eq!(round_tripped, original);
}
