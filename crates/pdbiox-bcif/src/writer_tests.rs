use super::*;
use crate::BinaryDocument;
use pdbiox_cif::{CifValue, DataBlock};
use pdbiox_core::io::Limits;
use pdbiox_core::span::ByteSpan;

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
            .push(value, pdbiox_cif::lexer::Quoting::Bare);
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
