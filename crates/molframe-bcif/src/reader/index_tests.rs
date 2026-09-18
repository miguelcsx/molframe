use super::BcifOffsetIndex;
use crate::container::{EncodedBlock, EncodedCategory, EncodedColumn, EncodedFile};
use crate::{EncodedData, Encoding, encode_integers};
use molframe_core::{ByteWindow, ExecutionContext, InputBuffer, SourceBytes};

#[test]
fn index_records_payload_offsets_without_copying_the_payloads() {
    let encoded = encode_integers(&[1, 2, 3, 5, 8]).expect("integer encoding");
    let expected = encoded.data.clone();
    let bytes = file_bytes(encoded);
    let context = ExecutionContext::default();
    let mut source = InputBuffer::from_bytes(bytes.clone());
    let index = BcifOffsetIndex::build(&mut source, 31, &context).expect("bounded index");
    assert_eq!(index.categories.len(), 1);
    let category = &index.categories[0];
    assert_eq!(category.name, "_atom_site");
    assert_eq!(category.rows, 5);
    assert_eq!(category.columns.len(), 1);
    let column = &category.columns[0];
    assert_eq!(column.name, "id");
    let start = usize::try_from(column.data.payload.start).expect("payload start");
    let end = usize::try_from(column.data.payload.end).expect("payload end");
    assert_eq!(&bytes[start..end], expected);
    assert!(column.data.encoding.start < column.data.encoding.end);
    let encoding_start = usize::try_from(column.data.encoding.start).expect("encoding start");
    let encoding_end = usize::try_from(column.data.encoding.end).expect("encoding end");
    let encoding: Vec<Encoding> = rmp_serde::from_slice(&bytes[encoding_start..encoding_end])
        .expect("indexed encoding metadata");
    assert_eq!(encoding, encoded_integer_encoding());
    assert!(column.mask.is_none());
    assert_eq!(context.live_batches(), 0);
}

#[test]
fn a_multi_gigabyte_binary_payload_is_skipped_without_reading_its_body() {
    let payload_length = u64::from(u32::MAX) - 32;
    let (prefix, suffix) = sparse_container(payload_length);
    let mut source = SparseSource::new(prefix, payload_length, suffix);
    let context = ExecutionContext::default();
    let index = BcifOffsetIndex::build(&mut source, 64, &context).expect("sparse bounded index");
    let column = &index.categories[0].columns[0];
    assert_eq!(
        column.data.payload.end - column.data.payload.start,
        payload_length
    );
    assert!(column.data.payload.start < u64::from(u32::MAX));
    assert!(column.data.payload.end > u64::from(u32::MAX));
    assert!(source.bytes_served < 1_024);
}

fn encoded_integer_encoding() -> Vec<Encoding> {
    encode_integers(&[1, 2, 3, 5, 8])
        .expect("integer encoding")
        .encoding
}

fn file_bytes(data: EncodedData) -> Vec<u8> {
    let file = EncodedFile {
        version: "0.3.0".to_owned(),
        encoder: "index-test".to_owned(),
        data_blocks: vec![EncodedBlock {
            header: "test".to_owned(),
            categories: vec![EncodedCategory {
                name: "_atom_site".to_owned(),
                row_count: 5,
                columns: vec![EncodedColumn {
                    name: "id".to_owned(),
                    data,
                    mask: None,
                }],
            }],
        }],
    };
    crate::messagepack::consume_to_vec_named(file).expect("serialize index fixture")
}

fn sparse_container(payload_length: u64) -> (Vec<u8>, Vec<u8>) {
    let length = u32::try_from(payload_length).expect("MessagePack bin32 length");
    let encoding = rmp_serde::to_vec_named(&vec![Encoding::ByteArray {
        r#type: crate::DataType::Uint8,
    }])
    .expect("encoding metadata");
    let mut prefix = Vec::new();
    prefix.extend_from_slice(&[
        0x81, 0xaa, b'd', b'a', b't', b'a', b'B', b'l', b'o', b'c', b'k', b's', 0x91, 0x81, 0xaa,
        b'c', b'a', b't', b'e', b'g', b'o', b'r', b'i', b'e', b's', 0x91, 0x83, 0xa4, b'n', b'a',
        b'm', b'e', 0xaa, b'_', b'a', b't', b'o', b'm', b'_', b's', b'i', b't', b'e', 0xa8, b'r',
        b'o', b'w', b'C', b'o', b'u', b'n', b't', 0x01, 0xa7, b'c', b'o', b'l', b'u', b'm', b'n',
        b's', 0x91, 0x82, 0xa4, b'n', b'a', b'm', b'e', 0xa2, b'i', b'd', 0xa4, b'd', b'a', b't',
        b'a', 0x82, 0xa8, b'e', b'n', b'c', b'o', b'd', b'i', b'n', b'g',
    ]);
    prefix.extend_from_slice(&encoding);
    prefix.extend_from_slice(&[0xa4, b'd', b'a', b't', b'a', 0xc6]);
    prefix.extend_from_slice(&length.to_be_bytes());
    (prefix, Vec::new())
}

struct SparseSource {
    prefix: Vec<u8>,
    payload_length: u64,
    suffix: Vec<u8>,
    buffer: Vec<u8>,
    bytes_served: u64,
}

impl SparseSource {
    fn new(prefix: Vec<u8>, payload_length: u64, suffix: Vec<u8>) -> Self {
        Self {
            prefix,
            payload_length,
            suffix,
            buffer: Vec::new(),
            bytes_served: 0,
        }
    }

    fn payload_start(&self) -> u64 {
        u64::try_from(self.prefix.len()).expect("prefix length")
    }
}

impl SourceBytes for SparseSource {
    fn window(
        &mut self,
        start: u64,
        length: usize,
    ) -> Result<ByteWindow<'_>, molframe_core::Diagnostic> {
        self.buffer.clear();
        let payload_start = self.payload_start();
        let payload_end = payload_start + self.payload_length;
        for offset in start..start.saturating_add(length as u64) {
            let value = if offset < payload_start {
                usize::try_from(offset)
                    .ok()
                    .and_then(|index| self.prefix.get(index).copied())
            } else if offset < payload_end {
                Some(0)
            } else {
                let local = usize::try_from(offset - payload_end).ok();
                local.and_then(|index| self.suffix.get(index).copied())
            };
            let Some(value) = value else { break };
            self.buffer.push(value);
        }
        self.bytes_served = self
            .bytes_served
            .saturating_add(u64::try_from(self.buffer.len()).expect("served bytes"));
        Ok(ByteWindow::new(start, &self.buffer))
    }

    fn len_hint(&self) -> Option<u64> {
        Some(self.payload_start() + self.payload_length + self.suffix.len() as u64)
    }
}
