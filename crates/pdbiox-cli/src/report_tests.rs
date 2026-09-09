use super::*;

#[test]
fn an_object_renders_its_fields_in_the_order_they_were_added() {
    let mut object = Json::new();
    object.text("id", "1ABC").number("atoms", 1_960);
    assert_eq!(object.finish(), r#"{"id":"1ABC","atoms":1960}"#);
}

#[test]
fn an_empty_object_renders_as_one() {
    assert_eq!(Json::new().finish(), "{}");
    assert_eq!(json_array(&[]), "[]");
}

#[test]
fn characters_a_string_cannot_carry_literally_are_escaped() {
    let mut object = Json::new();
    object.text("title", "a \"quoted\" \\ line\nbreak\ttab");
    assert_eq!(
        object.finish(),
        r#"{"title":"a \"quoted\" \\ line\nbreak\ttab"}"#
    );
}

#[test]
fn a_control_character_is_escaped_by_its_code_point() {
    let mut object = Json::new();
    object.text("odd", "\u{7}");
    assert_eq!(object.finish(), "{\"odd\":\"\\u0007\"}");
}

#[test]
fn a_list_of_objects_renders_as_an_array() {
    let mut first = Json::new();
    first.text("chain", "A");
    let mut second = Json::new();
    second.text("chain", "B");
    assert_eq!(
        json_array(&[first.finish(), second.finish()]),
        r#"[{"chain":"A"},{"chain":"B"}]"#
    );
}

#[test]
fn json_lines_keeps_one_object_per_line() {
    let context = Context {
        format: OutputKind::JsonLines,
        quiet: true,
        color: false,
        mode: pdbiox::ParseMode::Strict,
        policy: Box::leak(Box::new(pdbiox::AnalysisPolicy::default())),
        output: None,
        provenance: None,
        ccd: None,
        ccd_version: None,
        execution: Box::leak(Box::new(pdbiox::core::ExecutionContext::default())),
        missing_element_policy: pdbiox::MissingElementPolicy::PreserveUnknown,
        residue_boundary_policy: pdbiox::AmbiguousResidueBoundaryPolicy::Reject,
    };
    assert_eq!(
        context.json_records(&["{\"row\":1}".to_owned(), "{\"row\":2}".to_owned()]),
        "{\"row\":1}\n{\"row\":2}"
    );
}

#[test]
fn csv_quotes_delimiters_quotes_and_newlines_without_changing_row_order() {
    let mut table = Table::new(',', &["id", "title"]);
    table.row(["1ABC", "a, \"quoted\" title"]);
    table.row(["2DEF", "two\nlines"]);
    assert_eq!(
        table.finish(),
        "id,title\n1ABC,\"a, \"\"quoted\"\" title\"\n2DEF,\"two\nlines\""
    );
}

#[test]
fn json_results_embed_policy_and_read_provenance() {
    let context = Context {
        format: OutputKind::Json,
        quiet: true,
        color: false,
        mode: pdbiox::ParseMode::Strict,
        policy: Box::leak(Box::new(pdbiox::AnalysisPolicy::default())),
        output: None,
        provenance: None,
        ccd: None,
        ccd_version: None,
        execution: Box::leak(Box::new(pdbiox::core::ExecutionContext::default())),
        missing_element_policy: pdbiox::MissingElementPolicy::InferFromAtomName,
        residue_boundary_policy: pdbiox::AmbiguousResidueBoundaryPolicy::InferFromFileOrder,
    };
    let rendered = context.with_embedded_provenance("{\"atoms\":2}");
    assert!(rendered.contains("\"_provenance\""));
    assert!(rendered.contains("infer-from-atom-name"));
}

#[test]
fn tsv_only_quotes_fields_that_need_it() {
    let mut table = Table::new('\t', &["chain", "atoms"]);
    table.row(["A\tB", "42"]);
    assert_eq!(table.finish(), "chain\tatoms\n\"A\tB\"\t42");
}

#[test]
fn incremental_json_rows_are_published_as_one_valid_document() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path: &'static Path = Box::leak(Box::new(directory.path().join("rows.json")));
    let context = test_context(OutputKind::Json, Some(path));
    let mut rows = RowWriter::new(context, &["id", "title"]).expect("open rows");
    rows.row(["1", "first"]).expect("first row");
    rows.row(["2", "second"]).expect("second row");
    rows.finish().expect("finish rows");
    let output = std::fs::read_to_string(path).expect("published output");
    assert!(output.starts_with(
        "{\"result\":[{\"id\":\"1\",\"title\":\"first\"},{\"id\":\"2\",\"title\":\"second\"}]"
    ));
    assert!(output.contains("\"_provenance\""));
}

#[test]
fn an_unfinished_row_file_does_not_replace_the_destination() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path: &'static Path = Box::leak(Box::new(directory.path().join("rows.csv")));
    std::fs::write(path, "previous\n").expect("seed destination");
    let context = test_context(OutputKind::Csv, Some(path));
    let mut rows = RowWriter::new(context, &["id"]).expect("open rows");
    rows.row(["new"]).expect("row");
    drop(rows);
    assert_eq!(
        std::fs::read_to_string(path).expect("unchanged destination"),
        "previous\n"
    );
}

#[test]
fn invalid_row_width_writes_nothing_and_the_sink_remains_usable() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path: &'static Path = Box::leak(Box::new(directory.path().join("rows.csv")));
    let context = test_context(OutputKind::Csv, Some(path));
    let mut rows = RowWriter::new(context, &["id", "value"]).expect("open rows");
    assert!(rows.row(["too", "many", "columns"]).is_err());
    rows.row(["1", "kept"]).expect("valid row");
    rows.finish().expect("finish rows");
    assert_eq!(
        std::fs::read_to_string(path)
            .expect("published output")
            .lines()
            .take(2)
            .collect::<Vec<_>>(),
        ["id,value", "1,kept"]
    );
}

fn test_context(format: OutputKind, output: Option<&'static Path>) -> Context {
    Context {
        format,
        quiet: true,
        color: false,
        mode: pdbiox::ParseMode::Strict,
        policy: Box::leak(Box::new(pdbiox::AnalysisPolicy::default())),
        output,
        provenance: None,
        ccd: None,
        ccd_version: None,
        execution: Box::leak(Box::new(pdbiox::core::ExecutionContext::default())),
        missing_element_policy: pdbiox::MissingElementPolicy::PreserveUnknown,
        residue_boundary_policy: pdbiox::AmbiguousResidueBoundaryPolicy::Reject,
    }
}
