use super::*;

fn entry() -> Entry {
    Entry {
        name: "t".to_owned(),
        version: "1".to_owned(),
        url: "https://example.org/t".to_owned(),
        licence: "MIT".to_owned(),
        sha256: "AB".repeat(32),
        extension: "t".to_owned(),
        source: None,
        converter: None,
        append: Vec::new(),
    }
}

#[test]
fn a_written_table_reads_back_with_its_identity_digest_and_rows() {
    let text = table_text(&entry(), &["a\t1".to_owned(), "b\t2".to_owned()]);
    assert!(text.starts_with("# SPDX-License-Identifier: MIT\n"));
    let stored = read_back(&text).expect("valid");
    assert_eq!((stored.name.as_str(), stored.rows), ("t", 2));
    assert_eq!(stored.digest, Some("ab".repeat(32)));
}

#[test]
fn appended_files_contribute_rows_but_not_headers_or_comments() {
    let rows = bare_rows("# SPDX\nformat: 1\nname: x\n\nlink\t*\tC-N\t1.3\t0.01\n");
    assert_eq!(rows, ["link\t*\tC-N\t1.3\t0.01"]);
}
