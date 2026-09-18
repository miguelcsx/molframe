use super::*;

#[test]
fn a_bare_sentinel_is_a_sentinel_and_a_quoted_one_is_a_string() {
    assert_eq!(CifValue::parse(".", Quoting::Bare), CifValue::Inapplicable);
    assert_eq!(CifValue::parse("?", Quoting::Bare), CifValue::Unknown);
    assert_eq!(
        CifValue::parse(".", Quoting::Single),
        CifValue::Text(".".into())
    );
    assert_eq!(
        CifValue::parse("?", Quoting::Double),
        CifValue::Text("?".into())
    );
}

#[test]
fn a_whole_number_stays_whole_rather_than_becoming_a_float() {
    assert_eq!(CifValue::parse("42", Quoting::Bare), CifValue::Integer(42));
    assert_eq!(CifValue::parse("-7", Quoting::Bare), CifValue::Integer(-7));
    let large = 9_007_199_254_740_993_i64;
    assert_eq!(
        CifValue::parse(&large.to_string(), Quoting::Bare).as_integer(),
        Some(large),
        "a serial past the exact range of a float must survive"
    );
}

#[test]
fn a_number_with_a_fraction_is_read_as_one() {
    assert_eq!(CifValue::parse("1.5", Quoting::Bare), CifValue::Float(1.5));
    assert_eq!(CifValue::parse("42", Quoting::Bare).as_float(), Some(42.0));
}

#[test]
fn a_quoted_number_is_text_because_the_file_said_so() {
    assert_eq!(
        CifValue::parse("42", Quoting::Single),
        CifValue::Text("42".into())
    );
}

#[test]
fn an_identifier_reads_whether_it_was_written_as_a_word_or_a_number() {
    assert_eq!(
        CifValue::parse("A", Quoting::Bare)
            .as_identifier()
            .as_deref(),
        Some("A")
    );
    assert_eq!(
        CifValue::parse("2", Quoting::Bare)
            .as_identifier()
            .as_deref(),
        Some("2")
    );
    assert_eq!(CifValue::parse(".", Quoting::Bare).as_identifier(), None);
}

#[test]
fn a_sentinel_is_not_a_recorded_value() {
    assert!(CifValue::parse("x", Quoting::Bare).is_recorded());
    assert!(!CifValue::parse(".", Quoting::Bare).is_recorded());
    assert!(!CifValue::parse("?", Quoting::Bare).is_recorded());
}

#[test]
fn categories_and_items_keep_the_order_the_file_gave_them() {
    let mut block = DataBlock::new("test");
    let category = block.category_mut("atom_site", ByteSpan::default());
    for item in ["group_PDB", "id", "type_symbol"] {
        category
            .column_mut(item)
            .push(CifValue::Integer(1), Quoting::Bare);
    }
    let items: Vec<_> = category.items().collect();
    assert_eq!(items, ["group_PDB", "id", "type_symbol"]);
    assert_eq!(category.row_count(), 1);
}

#[test]
fn a_category_the_library_has_no_interpretation_for_is_kept_all_the_same() {
    let mut block = DataBlock::new("test");
    block
        .category_mut("my_lab_custom_category", ByteSpan::default())
        .column_mut("value")
        .push(CifValue::Text("kept".into()), Quoting::Bare);
    assert_eq!(
        block.category("my_lab_custom_category").map(Category::len),
        Some(1)
    );
    assert_eq!(
        block
            .category("my_lab_custom_category")
            .and_then(|c| c.text("value", 0)),
        Some("kept")
    );
}

#[test]
fn declined_shapes_still_parse_through_the_general_path() {
    let value = CifValue::parse("1.5e-3", Quoting::Bare);
    assert_eq!(value.as_float(), Some(0.0015));
    let value = CifValue::parse("1e5", Quoting::Bare);
    assert_eq!(value.as_float(), Some(100_000.0));
}
