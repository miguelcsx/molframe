//! A read can record the SHA-256 of the bytes it parsed, and only when asked.

use crate::{ReadOptions, read_bytes};
use molframe_core::contract::ContentDigest;

const PDB: &str =
    "ATOM      1  N   GLY A   1       1.000   1.000   1.000  1.00 10.00           N\nEND\n";

fn read(options: &ReadOptions) -> crate::Structure {
    match read_bytes(PDB.as_bytes().to_vec(), Some("t.pdb"), options) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("the fixture reads: {findings:?}"),
    }
}

#[test]
fn a_read_that_asks_for_the_digest_records_the_digest_of_its_bytes() {
    let structure = read(&ReadOptions::new().digest_input(true));
    let entry = &structure.metadata();
    assert_eq!(entry.input_sha256, Some(ContentDigest::of(PDB.as_bytes())));
    assert_eq!(entry.input_name.as_deref(), Some("t.pdb"));
}

#[test]
fn a_read_that_does_not_ask_records_the_name_and_no_digest() {
    let structure = read(&ReadOptions::new());
    let entry = &structure.metadata();
    assert_eq!(entry.input_sha256, None);
    assert_eq!(entry.input_name.as_deref(), Some("t.pdb"));
}
