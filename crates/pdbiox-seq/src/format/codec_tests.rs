use super::*;

#[test]
fn every_record_format_round_trips_through_common_dispatch() {
    let cases = [
        (SequenceFormat::Fasta, ">a\nAC-G\n>b\nACTG\n"),
        (
            SequenceFormat::Stockholm,
            "# STOCKHOLM 1.0\na AC-G\nb ACTG\n//\n",
        ),
        (SequenceFormat::Clustal, "CLUSTAL W\n\na AC-G\nb ACTG\n"),
        (SequenceFormat::Phylip, "2 4\na AC-G\nb ACTG\n"),
        (SequenceFormat::A2m, ">a\nAC-gT\n>b\nACT-\n"),
        (SequenceFormat::A3m, ">a\nACg-\n>b\nAC-\n"),
    ];
    for (format, source) in cases {
        let Ok(document) = read_sequence(source, format) else {
            panic!("{format:?} should read");
        };
        let Ok(rendered) = write_sequence(&document, format) else {
            panic!("{format:?} should write");
        };
        assert!(read_sequence(&rendered, format).is_ok(), "{format:?}");
    }
}

#[test]
fn fastq_and_newick_keep_their_distinct_document_kinds() {
    let fastq = read_sequence("@x\nAC\n+\n!!\n", SequenceFormat::Fastq);
    assert!(matches!(fastq, Ok(SequenceDocument::Fastq(_))));
    let tree = read_sequence("(a:1,b:1);", SequenceFormat::Newick);
    assert!(matches!(tree, Ok(SequenceDocument::Tree(_))));
}

#[test]
fn malformed_dimensions_and_wrong_document_variants_are_refused() {
    assert!(read_sequence("2 5\na AC\nb AC\n", SequenceFormat::Phylip).is_err());
    let tree = SequenceDocument::Tree(Tree::Leaf { name: "x".into() });
    assert!(write_sequence(&tree, SequenceFormat::Fasta).is_err());
}
