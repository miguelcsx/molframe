use super::*;

fn buffer(text: &str) -> InputBuffer {
    InputBuffer::from_bytes(text.as_bytes().to_vec())
}

#[test]
fn structured_text_is_recognised_by_its_block_header() {
    assert!(Format::Mmcif.recognises(b"data_1ABC\n_entry.id 1ABC\n"));
    assert!(Format::Mmcif.recognises(b"# a comment\ndata_1ABC\n"));
    assert!(!Format::Mmcif.recognises(b"ATOM      1  N   ALA A   1\n"));
}

#[test]
fn fixed_column_files_are_recognised_by_their_record_names() {
    assert!(Format::Pdb.recognises(b"ATOM      1  N   ALA A   1\n"));
    assert!(Format::Pdb.recognises(b"HEADER    HYDROLASE\nATOM      1\n"));
    assert!(Format::Pdb.recognises(b"HETATM    1 ZN    ZN A 100\n"));
    assert!(!Format::Pdb.recognises(b"data_1ABC\n"));
}

#[test]
fn variant_names_and_autodock_markers_select_their_formats() {
    assert_eq!(Format::from_name("ligand.pqr"), Some(Format::Pqr));
    assert_eq!(Format::from_name("ligand.pdbqt.gz"), Some(Format::Pdbqt));
    assert!(Format::Pdbqt.recognises(b"ROOT\nTORSDOF 0\n"));
    assert_eq!(Format::Pqr.name(), "pqr");
}

#[test]
fn content_decides_before_the_file_name_does() {
    let mislabelled = buffer("data_1ABC\n_entry.id 1ABC\n");
    let detected = Format::detect(Format::Auto, &mislabelled, Some("entry.pdb"));
    assert_eq!(detected.ok(), Some(Format::Mmcif));
}

#[test]
fn the_file_name_decides_when_the_content_recognises_nothing() {
    let opaque = buffer("something else entirely\n");
    assert_eq!(
        Format::detect(Format::Auto, &opaque, Some("x.cif")).ok(),
        Some(Format::Mmcif)
    );
    assert_eq!(
        Format::detect(Format::Auto, &opaque, Some("x.pdb")).ok(),
        Some(Format::Pdb)
    );
}

#[test]
fn a_compression_suffix_does_not_hide_the_format_suffix_under_it() {
    assert_eq!(Format::from_name("1abc.cif.gz"), Some(Format::Mmcif));
    assert_eq!(Format::from_name("1abc.pdb.zst"), Some(Format::Pdb));
    assert_eq!(Format::from_name("1abc.CIF"), Some(Format::Mmcif));
    assert_eq!(Format::from_name("1abc.txt"), None);
}

#[test]
fn an_unrecognisable_input_names_what_was_tried_rather_than_guessing() {
    let opaque = buffer("nothing recognisable\n");
    let refused = Format::detect(Format::Auto, &opaque, Some("x.txt"));
    assert_eq!(
        refused.err().map(|finding| finding.code()),
        Some(Code::E1001)
    );
}

#[test]
fn a_named_format_is_taken_at_its_word_without_sniffing() {
    let opaque = buffer("nothing recognisable\n");
    assert_eq!(
        Format::detect(Format::Pdb, &opaque, Some("x.cif")).ok(),
        Some(Format::Pdb)
    );
}

#[test]
fn parse_mode_maps_onto_how_severe_a_finding_must_be_to_stop_a_read() {
    assert_eq!(ParseMode::Strict.strictness(), Strictness::Strict);
    assert_eq!(ParseMode::Permissive.strictness(), Strictness::Medium);
    assert_eq!(ParseMode::Recover.strictness(), Strictness::Loose);
    assert_eq!(ParseMode::default(), ParseMode::Permissive);
}

#[test]
fn a_filter_that_overrides_nothing_writes_everything() {
    use crate::index::AtomIndex;
    assert!(SelectAll.accept_atom(AtomIndex::new(0)));
    assert!(SelectAll.accept_model(crate::index::ModelIndex::new(9)));
}

#[test]
fn mmtf_messagepack_is_distinguished_from_binary_cif() {
    let bytes = b"\x81\xabmmtfVersion\xa51.0.0";
    assert!(Format::Mmtf.recognises(bytes));
    assert!(!Format::BinaryCif.recognises(bytes));
    assert_eq!(Format::from_name("entry.mmtf.gz"), Some(Format::Mmtf));
}

#[test]
fn pdbml_is_recognized_by_xml_root_and_suffix() {
    let bytes = br#"<?xml version="1.0"?><PDBx:datablock datablockName="X"/>"#;
    assert!(Format::Pdbml.recognises(bytes));
    assert_eq!(Format::from_name("entry.xml.gz"), Some(Format::Pdbml));
}
