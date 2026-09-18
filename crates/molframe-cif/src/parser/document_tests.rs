use super::*;
use std::fmt::Write as _;
use std::sync::Arc;

fn document(text: &str) -> Document {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    match parse(&input) {
        Ok((document, _)) => document,
        Err(findings) => panic!("parse failed: {findings:?}"),
    }
}

fn findings(text: &str) -> Vec<Code> {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let reported = match parse(&input) {
        Ok((_, findings)) | Err(findings) => findings,
    };
    reported.iter().map(Diagnostic::code).collect()
}

#[test]
fn a_tag_and_value_pair_becomes_one_row_of_one_category() {
    let document = document("data_x\n_entry.id 1ABC\n");
    let id = document
        .first_block()
        .and_then(|block| block.category("entry"))
        .and_then(|entry| entry.text("id", 0));
    assert_eq!(id, Some("1ABC"));
}

#[test]
fn a_loop_becomes_as_many_rows_as_it_has_value_groups() {
    let document =
        document("data_x\nloop_\n_atom_site.id\n_atom_site.type_symbol\n1 N\n2 C\n3 O\n");
    let category = document
        .first_block()
        .and_then(|block| block.category("atom_site"));
    assert_eq!(category.map(Category::row_count), Some(3));
    assert_eq!(category.and_then(|c| c.text("type_symbol", 2)), Some("O"));
    assert_eq!(
        category
            .and_then(|c| c.value("id", 1))
            .and_then(CifValue::as_integer),
        Some(2)
    );
}

#[test]
fn repeated_text_values_share_one_allocation() {
    let document = document("data_x\nloop_\n_a.name\nALA\nALA\n'ALA'\n");
    let Some(category) = document.first_block().and_then(|block| block.category("a")) else {
        panic!("expected a category")
    };
    let Some(CifValue::Text(first)) = category.value("name", 0) else {
        panic!("expected text")
    };
    let Some(CifValue::Text(second)) = category.value("name", 1) else {
        panic!("expected text")
    };
    let Some(CifValue::Text(quoted)) = category.value("name", 2) else {
        panic!("expected text")
    };
    assert!(Arc::ptr_eq(first, second));
    assert!(Arc::ptr_eq(first, quoted));
}

#[test]
fn an_item_name_splits_into_its_category_and_its_item() {
    assert_eq!(split_tag("_atom_site.id"), ("atom_site", "id"));
    assert_eq!(split_tag("atom_site.id"), ("atom_site", "id"));
    assert_eq!(split_tag("_legacy"), ("legacy", "legacy"));
}

#[test]
fn a_loop_row_that_is_short_is_reported_with_what_it_should_have_had() {
    let codes = findings("data_x\nloop_\n_a.one\n_a.two\n_a.three\n1 2 3\n4 5\n");
    assert!(codes.contains(&Code::E1103), "got {codes:?}");
}

#[test]
fn a_file_with_no_block_header_is_refused() {
    let codes = findings("_entry.id 1ABC\n");
    assert!(codes.contains(&Code::E1106), "got {codes:?}");
}

#[test]
fn several_blocks_are_all_kept() {
    let document = document("data_one\n_a.b 1\ndata_two\n_a.b 2\n");
    assert_eq!(document.len(), 2);
    let names: Vec<_> = document.blocks().map(DataBlock::name).collect();
    assert_eq!(names, ["one", "two"]);
}

#[test]
fn rows_walk_a_category_one_row_at_a_time() {
    let document = document("data_x\nloop_\n_a.id\n_a.name\n1 first\n2 second\n");
    let Some(category) = document.first_block().and_then(|b| b.category("a")) else {
        panic!("expected a category")
    };
    let mut rows = Rows::new(category);
    assert_eq!(rows.text("name"), Some("first"));
    assert_eq!(rows.integer("id"), Some(1));
    assert!(rows.advance());
    assert_eq!(rows.text("name"), Some("second"));
    assert!(!rows.advance());
}

#[test]
fn an_item_a_row_does_not_carry_reads_as_absent_rather_than_as_a_guess() {
    let document = document("data_x\n_a.present yes\n");
    let Some(category) = document.first_block().and_then(|b| b.category("a")) else {
        panic!("expected a category")
    };
    let rows = Rows::new(category);
    assert_eq!(rows.text("absent"), None);
    assert!(!rows.is_recorded("absent"));
}

#[derive(Clone, Debug, PartialEq)]
enum OwnedScalar {
    Inapplicable,
    Unknown,
    Text(String),
    Integer(i64),
    Float(f64),
}

#[derive(Default)]
struct EventProbe {
    values: Vec<(String, String, OwnedScalar)>,
    rows: usize,
}

impl CifEventSink for EventProbe {
    type Output = Self;

    fn block(&mut self, _name: &str) {}

    fn value(&mut self, category: &str, item: &str, value: CifScalar<'_>, _span: ByteSpan) {
        let owned = match value {
            CifScalar::Inapplicable => OwnedScalar::Inapplicable,
            CifScalar::Unknown => OwnedScalar::Unknown,
            CifScalar::Text(value) => OwnedScalar::Text(value.to_owned()),
            CifScalar::Integer(value) => OwnedScalar::Integer(value),
            CifScalar::Float(value) => OwnedScalar::Float(value),
        };
        self.values
            .push((category.to_owned(), item.to_owned(), owned));
    }

    fn end_row(&mut self) {
        self.rows += 1;
    }

    fn finish(self) -> Self::Output {
        self
    }
}

#[test]
fn event_projection_matches_scalar_loop_sentinel_and_multiline_semantics() {
    let source = "data_x\n_entry.id 'model one'\nloop_\n_a.id\n_a.score\n_a.note\n\
1 2.5 .\n2 ?\n;line one\nline two\n;\n";
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    let (_, document_findings) = parse(&input).expect("document oracle parses");
    let (events, event_findings) =
        parse_events(&input, EventProbe::default()).expect("event projection parses");
    assert_eq!(event_findings, document_findings);
    assert_eq!(events.rows, 3);
    assert_eq!(
        events.values,
        vec![
            (
                "entry".to_owned(),
                "id".to_owned(),
                OwnedScalar::Text("model one".to_owned())
            ),
            ("a".to_owned(), "id".to_owned(), OwnedScalar::Integer(1)),
            ("a".to_owned(), "score".to_owned(), OwnedScalar::Float(2.5)),
            ("a".to_owned(), "note".to_owned(), OwnedScalar::Inapplicable),
            ("a".to_owned(), "id".to_owned(), OwnedScalar::Integer(2)),
            ("a".to_owned(), "score".to_owned(), OwnedScalar::Unknown),
            (
                "a".to_owned(),
                "note".to_owned(),
                OwnedScalar::Text("line one\nline two".to_owned())
            ),
        ]
    );
}

#[test]
fn event_projection_and_document_report_the_same_malformed_loop() {
    let source = "data_x\nloop_\n_a.one\n_a.two\n1 2\n3\n";
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    let document = match parse(&input) {
        Ok((_, findings)) | Err(findings) => findings,
    };
    let events = match parse_events(&input, EventProbe::default()) {
        Ok((_, findings)) | Err(findings) => findings,
    };
    assert_eq!(
        events.iter().map(Diagnostic::code).collect::<Vec<_>>(),
        document.iter().map(Diagnostic::code).collect::<Vec<_>>()
    );
}

#[test]
fn event_projection_scales_without_materialising_cells() {
    let mut source = String::from("data_x\nloop_\n_a.id\n_a.value\n");
    for row in 0..10_000 {
        let _ = writeln!(source, "{row} repeated");
    }
    let input = InputBuffer::from_bytes(source.into_bytes());
    let (events, findings) =
        parse_events(&input, EventCounter::default()).expect("large event loop parses");
    assert!(findings.is_empty());
    assert_eq!(events.rows, 10_000);
    assert_eq!(events.values, 20_000);
}

#[derive(Default)]
struct EventCounter {
    values: usize,
    rows: usize,
}

impl CifEventSink for EventCounter {
    type Output = Self;

    fn block(&mut self, _name: &str) {}

    fn value(&mut self, _category: &str, _item: &str, _value: CifScalar<'_>, _span: ByteSpan) {
        self.values += 1;
    }

    fn end_row(&mut self) {
        self.rows += 1;
    }

    fn finish(self) -> Self::Output {
        self
    }
}
