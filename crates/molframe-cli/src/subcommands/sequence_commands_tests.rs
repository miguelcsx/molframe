use super::*;

#[test]
fn format_conversion_delegates_to_native_parsers_and_writers() {
    let document = molframe::sequence::read_sequence(
        ">a\nACD\n>b\nA-D\n",
        molframe::sequence::SequenceFormat::Fasta,
    )
    .unwrap_or_else(|error| panic!("parse failed: {error}"));
    let rendered =
        molframe::sequence::write_sequence(&document, molframe::sequence::SequenceFormat::Clustal)
            .unwrap_or_else(|error| panic!("write failed: {error}"));
    assert!(rendered.starts_with("CLUSTAL"));
    assert!(rendered.contains('a'));
    assert!(rendered.contains('b'));
}

#[test]
fn fastq_output_refuses_to_invent_quality_scores() {
    let document =
        molframe::sequence::read_sequence(">a\nACD\n", molframe::sequence::SequenceFormat::Fasta)
            .unwrap_or_else(|error| panic!("parse failed: {error}"));
    let error =
        molframe::sequence::write_sequence(&document, molframe::sequence::SequenceFormat::Fastq)
            .expect_err("FASTA has no quality scores");
    assert!(error.to_string().contains("incompatible"));
}

#[test]
fn distance_matrix_requires_a_square_body() {
    let error = parse_distance_matrix("a\tb\n0\t1\n").expect_err("one row is not square");
    assert!(error.contains("square"));
}

#[test]
fn alignment_projection_preserves_native_correspondence() {
    let alignment =
        match molframe::sequence::global(b"AC", b"A", molframe::sequence::Scoring::simple()) {
            Ok(alignment) => alignment,
            Err(error) => panic!("alignment failed: {error}"),
        };
    let (left, right) = aligned_sequences(b"AC", b"A", &alignment)
        .unwrap_or_else(|error| panic!("projection failed: {error}"));
    assert_eq!(left.len(), right.len());
    assert_eq!(left, b"AC");
    assert_eq!(right, b"A-");
}
