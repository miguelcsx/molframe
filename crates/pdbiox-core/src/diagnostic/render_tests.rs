use super::*;
use crate::diagnostic::Code;
use crate::span::{ByteSpan, Position};

#[test]
fn a_finding_renders_its_severity_word_code_and_message() {
    let finding = Diagnostic::new(Code::E1106);
    let text = Rendered::new(&finding).to_string();
    assert!(text.starts_with("error[PDBIOX-E1106]: missing data_ block header"));
    assert!(text.contains("help: add a data_ line"));
}

#[test]
fn an_informational_finding_renders_as_a_note_not_an_error() {
    let finding = Diagnostic::new(Code::W2001);
    assert!(
        Rendered::new(&finding)
            .to_string()
            .starts_with("note[PDBIOX-W2001]")
    );
}

#[test]
fn a_span_with_source_quotes_the_offending_line_and_underlines_it() {
    let source = b"data_test\nloop_\n_atom_site.id\n";
    let span = ByteSpan::new(Position::new(10, 2, 1), 15);
    let finding = Diagnostic::new(Code::E1103)
        .at(span)
        .in_category("atom_site");
    let text = Rendered::new(&finding)
        .with_source(source)
        .with_origin("t.cif")
        .to_string();

    assert!(text.contains("t.cif:2:1"));
    assert!(text.contains("loop_"));
    assert!(text.contains("^^^^^"));
    assert!(text.contains("category: atom_site"));
}

#[test]
fn rendering_without_colour_emits_no_escape_sequences() {
    let finding = Diagnostic::new(Code::E1103).at(ByteSpan::empty(Position::START));
    let text = Rendered::new(&finding).with_source(b"loop_\n").to_string();
    assert!(!text.contains('\u{1b}'));
}

#[test]
fn rendering_with_colour_emits_escape_sequences() {
    let finding = Diagnostic::new(Code::E1103);
    let text = Rendered::new(&finding).with_color(true).to_string();
    assert!(text.contains('\u{1b}'));
}

#[test]
fn a_span_beyond_the_source_is_skipped_rather_than_panicking() {
    let finding = Diagnostic::new(Code::E1103).at(ByteSpan::new(Position::new(500, 9, 1), 510));
    let text = Rendered::new(&finding).with_source(b"short\n").to_string();
    assert!(text.contains("PDBIOX-E1103"));
}

#[test]
fn line_lookup_finds_the_line_a_byte_offset_falls_on() {
    let source = b"one\ntwo\nthree";
    let line = line_containing(source, 5).map(|line| line.text);
    assert_eq!(line, Some(&b"two"[..]));
    assert_eq!(
        line_containing(source, 0).map(|line| line.text),
        Some(&b"one"[..])
    );
    assert_eq!(line_containing(source, 999).map(|line| line.text), None);
}
