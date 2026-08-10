use super::*;

#[test]
fn format_conversion_delegates_to_native_parsers_and_writers() {
    let document =
        pdbiox::seq::read_sequence(">a\nACD\n>b\nA-D\n", pdbiox::seq::SequenceFormat::Fasta)
            .unwrap_or_else(|error| panic!("parse failed: {error}"));
    let rendered = pdbiox::seq::write_sequence(&document, pdbiox::seq::SequenceFormat::Clustal)
        .unwrap_or_else(|error| panic!("write failed: {error}"));
    assert!(rendered.starts_with("CLUSTAL"));
    assert!(rendered.contains('a'));
    assert!(rendered.contains('b'));
}

#[test]
fn fastq_output_refuses_to_invent_quality_scores() {
    let document = pdbiox::seq::read_sequence(">a\nACD\n", pdbiox::seq::SequenceFormat::Fasta)
        .unwrap_or_else(|error| panic!("parse failed: {error}"));
    let error = pdbiox::seq::write_sequence(&document, pdbiox::seq::SequenceFormat::Fastq)
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
    let alignment = match pdbiox::seq::global(b"AC", b"A", pdbiox::seq::Scoring::simple()) {
        Ok(alignment) => alignment,
        Err(error) => panic!("alignment failed: {error}"),
    };
    let (left, right) = aligned_sequences(b"AC", b"A", &alignment)
        .unwrap_or_else(|error| panic!("projection failed: {error}"));
    assert_eq!(left.len(), right.len());
    assert_eq!(left, b"AC");
    assert_eq!(right, b"A-");
}
