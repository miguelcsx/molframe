use super::*;

#[cfg(feature = "pdb")]
use crate::PdbOptions;

const DIPEPTIDE: &str = "\
ATOM      1  N   GLY A   1      27.340  24.430   2.614  1.00 10.00           N
ATOM      2  CA  GLY A   1      26.266  25.413   2.842  1.00 11.00           C
END
";

#[test]
fn bytes_are_dispatched_to_the_reader_their_content_names() {
    let read = read_bytes(DIPEPTIDE.as_bytes().to_vec(), None, &ReadOptions::new());
    assert_eq!(
        read.map(|(structure, _)| structure.atom_count()).ok(),
        Some(2)
    );
}

#[cfg(feature = "mmcif")]
#[test]
fn structured_text_is_dispatched_to_the_reader_its_content_names() {
    let text = b"data_1ABC\n_entry.id 1ABC\n".to_vec();
    // No coordinates, so the read is refused — but by the reader that owns the
    // format, which is the dispatch this checks.
    let refused = read_bytes(text, None, &ReadOptions::new());
    let findings = refused.err().unwrap_or_default();
    assert_eq!(findings.first().map(Diagnostic::code), Some(Code::E2001));
}

#[test]
fn an_unreadable_path_reports_a_finding_rather_than_panicking() {
    assert!(read("no/such/file.pdb").is_err());
}

#[cfg(feature = "pdb")]
#[test]
fn a_structure_read_through_the_facade_writes_back_through_it() {
    let Ok((structure, _)) = read_bytes(DIPEPTIDE.as_bytes().to_vec(), None, &ReadOptions::new())
    else {
        panic!("expected a structure")
    };
    let written = write_pdb(&structure, &PdbOptions::new());
    assert!(written.unwrap_or_default().contains("ATOM"));
}
