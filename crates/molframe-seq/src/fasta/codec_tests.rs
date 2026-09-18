use super::{FastaRecord, parse_fasta, write_fasta};

#[test]
fn records_are_read_in_order_with_sequences_joined() {
    let text = ">a a description\nACGT\nGGCC\n>b\nTTTT\n";
    let records = parse_fasta(text);
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].id, "a");
    assert_eq!(records[0].description, "a description");
    assert_eq!(records[0].sequence, b"ACGTGGCC");
    assert_eq!(records[1].id, "b");
    assert!(records[1].description.is_empty());
    assert_eq!(records[1].sequence, b"TTTT");
}

#[test]
fn text_before_the_first_header_is_ignored() {
    let records = parse_fasta("junk\n>a\nACGT\n");
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].sequence, b"ACGT");
}

#[test]
fn writing_then_reading_returns_the_same_records() {
    let records = vec![
        FastaRecord {
            id: "long".to_string(),
            description: "wrapped".to_string(),
            // Long enough to force wrapping across several lines.
            sequence: b"ACGT".iter().cycle().take(150).copied().collect(),
        },
        FastaRecord {
            id: "short".to_string(),
            description: String::new(),
            sequence: b"MKV".to_vec(),
        },
    ];
    let round_tripped = parse_fasta(&write_fasta(&records));
    assert_eq!(round_tripped, records);
}
