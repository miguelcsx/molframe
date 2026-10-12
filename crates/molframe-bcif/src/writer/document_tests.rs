use super::*;
use crate::BinaryDocument;
use molframe_cif::{CifValue, DataBlock};
use molframe_core::io::Limits;
use molframe_core::span::ByteSpan;
use std::sync::Arc;

#[test]
fn document_round_trip_keeps_columns_values_masks_and_order() {
    let mut source = Document::new();
    let mut block = DataBlock::new("test");
    let category = block.category_mut("numbers", ByteSpan::default());
    for value in [
        CifValue::Integer(7),
        CifValue::Unknown,
        CifValue::Integer(9),
    ] {
        category
            .column_mut("value")
            .push(value, molframe_cif::lexer::Quoting::Bare);
    }
    source.push(block);

    let bytes = write_document(&source).expect("document writes");
    let decoded = BinaryDocument::parse(&bytes, Limits::default())
        .expect("container opens")
        .to_document()
        .expect("document decodes");
    let values = decoded
        .first_block()
        .and_then(|block| block.category("numbers"))
        .and_then(|category| category.column("value"))
        .expect("column exists");
    assert_eq!(values.get(0), Some(&CifValue::Integer(7)));
    assert_eq!(values.get(1), Some(&CifValue::Unknown));
    assert_eq!(values.get(2), Some(&CifValue::Integer(9)));
}

#[test]
fn decoded_document_shares_repeated_text_allocations() {
    let mut source = Document::new();
    let mut block = DataBlock::new("test");
    let column = block
        .category_mut("labels", ByteSpan::default())
        .column_mut("value");
    for _ in 0..3 {
        column.push(
            CifValue::Text(Arc::from("repeated")),
            molframe_cif::lexer::Quoting::Bare,
        );
    }
    source.push(block);

    let bytes = write_document(&source).expect("document writes");
    let decoded = BinaryDocument::parse(&bytes, Limits::default())
        .expect("container opens")
        .to_document()
        .expect("document decodes");
    let values = decoded
        .first_block()
        .and_then(|block| block.category("labels"))
        .and_then(|category| category.column("value"))
        .expect("column exists");
    let Some(CifValue::Text(first)) = values.get(0) else {
        panic!("first text exists")
    };
    let Some(CifValue::Text(second)) = values.get(1) else {
        panic!("second text exists")
    };
    assert!(Arc::ptr_eq(first, second));
}

#[test]
fn column_mixing_numbers_and_text_is_written_as_strings() {
    let mut source = Document::new();
    let mut block = DataBlock::new("test");
    let category = block.category_mut("pdbx_struct_assembly_gen", ByteSpan::default());
    for value in [
        CifValue::Integer(1),
        CifValue::Text(Arc::from("(1-60)")),
        CifValue::Float(2.5),
        CifValue::Unknown,
    ] {
        category
            .column_mut("oper_expression")
            .push(value, molframe_cif::lexer::Quoting::Bare);
    }
    source.push(block);

    let bytes = write_document(&source).expect("mixed column writes");
    let decoded = BinaryDocument::parse(&bytes, Limits::default())
        .expect("container opens")
        .to_document()
        .expect("document decodes");
    let column = decoded
        .first_block()
        .and_then(|block| block.category("pdbx_struct_assembly_gen"))
        .and_then(|category| category.column("oper_expression"))
        .expect("column exists");
    let text = |row| match column.get(row) {
        Some(CifValue::Text(text)) => Some(text.to_string()),
        _ => None,
    };
    assert_eq!(text(0).as_deref(), Some("1"));
    assert_eq!(text(1).as_deref(), Some("(1-60)"));
    assert_eq!(text(2).as_deref(), Some("2.5"));
    assert_eq!(column.get(3), Some(&CifValue::Unknown));
}
