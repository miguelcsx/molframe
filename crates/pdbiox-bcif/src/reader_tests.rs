use super::*;
use crate::write_document;
use pdbiox_cif::{CifValue, DataBlock};
use pdbiox_core::span::ByteSpan;

fn document() -> Document {
    let mut document = Document::new();
    let mut block = DataBlock::new("test");
    let category = block.category_mut("entry", ByteSpan::default());
    category.column_mut("id").push(
        CifValue::Text("1ABC".into()),
        pdbiox_cif::lexer::Quoting::Bare,
    );
    document.push(block);
    document
}

#[test]
fn opening_a_document_keeps_columns_lazy_until_requested() {
    let bytes = write_document(&document()).expect("document writes");
    let input = InputBuffer::from_bytes(bytes);
    let binary = read_document(&input, Limits::default()).expect("document opens");
    assert_eq!(binary.block_count(), 1);
    let entry = binary
        .category(0, "entry")
        .expect("category decodes")
        .expect("entry exists");
    assert_eq!(entry.text("id", 0), Some("1ABC"));
}
