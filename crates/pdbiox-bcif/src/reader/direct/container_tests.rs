use super::*;
use pdbiox_cif::{CifValue, DataBlock, Document};

#[test]
fn direct_container_payloads_borrow_the_input_buffer() {
    let mut document = Document::new();
    let mut block = DataBlock::new("borrowed");
    block
        .category_mut("atom_site", ByteSpan::default())
        .column_mut("type_symbol")
        .push(CifValue::Text("C".into()), Quoting::Bare);
    document.push(block);
    let bytes = crate::write_document(&document).expect("fixture encodes");
    let file = parse(&bytes, Limits::default()).expect("container parses");
    let payload = file.data_blocks[0].categories[0].columns[0].data.data;
    let input = bytes.as_ptr_range();
    let payload_start = payload.as_ptr();
    let payload_end = payload[payload.len()..].as_ptr();

    assert!(payload_start >= input.start);
    assert!(payload_end <= input.end);
}
