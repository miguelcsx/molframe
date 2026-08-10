use super::*;

#[test]
fn exact_and_bucketed_tables_return_identical_hits() {
    let sequences: [&[u8]; 2] = [b"ACAC", b"TTAC"];
    let options = |storage| KmerTableOptions {
        k: 2,
        pattern: None,
        storage,
    };
    let Ok(exact) = KmerTable::build(&sequences, &options(KmerStorage::Exact)) else {
        panic!("valid exact table");
    };
    let Ok(bucketed) = KmerTable::build(&sequences, &options(KmerStorage::Bucketed { buckets: 3 }))
    else {
        panic!("valid bucketed table");
    };
    assert_eq!(exact.query(b"AC"), bucketed.query(b"AC"));
    assert_eq!(exact.query(b"AC").len(), 3);
}

#[test]
fn a_spaced_table_reports_source_window_positions() {
    let Ok(pattern) = SeedPattern::new(&[true, false, true]) else {
        panic!("valid pattern");
    };
    let options = KmerTableOptions {
        k: 0,
        pattern: Some(pattern),
        storage: KmerStorage::Exact,
    };
    let Ok(table) = KmerTable::build(&[b"ABCDE"], &options) else {
        panic!("valid spaced table");
    };
    assert_eq!(
        table.query(b"AC"),
        &[KmerHit {
            sequence: 0,
            position: 0
        }]
    );
}
