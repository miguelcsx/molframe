use super::*;
use molframe_core::io::Format;

const ONE_ATOM: &str =
    "ATOM      1  CA  ALA A   1      11.104   6.134  -6.504  1.00  0.00           C\nEND\n";

fn one_atom() -> Structure {
    match read_bytes(
        ONE_ATOM.as_bytes().to_vec(),
        Some("a.pdb"),
        &ReadOptions::new(),
    ) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("read failed: {findings:?}"),
    }
}

#[test]
fn an_explicit_format_writes_to_a_name_that_selects_none() {
    let structure = one_atom();
    let path = std::env::temp_dir().join(format!("molframe-write-as-{}.dat", std::process::id()));
    let options = WriteOptions::canonical();
    let refused = write_with_options(&path, &structure, &options);
    assert!(refused.is_err(), "the name alone selects no writer");
    if let Err(findings) = write_as(&path, &structure, Format::Pdb, &options) {
        panic!("explicit write failed: {findings:?}")
    }
    let written = std::fs::read_to_string(&path);
    let _removed = std::fs::remove_file(&path);
    assert!(written.is_ok_and(|text| text.contains("ATOM")));
}

#[test]
fn a_bare_extension_names_its_format_without_regard_to_case() {
    assert_eq!(Format::from_extension("PDB"), Some(Format::Pdb));
    assert_eq!(Format::from_extension("bcif"), Some(Format::BinaryCif));
    assert_eq!(Format::from_extension("nope"), None);
}
