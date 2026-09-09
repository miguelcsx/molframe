use super::*;
use crate::codec::{DataType, Encoding};

#[test]
fn consuming_serialization_keeps_exact_named_bytes() {
    let file = EncodedFile {
        version: "0.3.0".to_owned(),
        encoder: "pdbiox test".to_owned(),
        data_blocks: vec![EncodedBlock {
            header: "TEST".to_owned(),
            categories: vec![EncodedCategory {
                name: "_atom_site".to_owned(),
                row_count: 3,
                columns: vec![
                    column("id", vec![1, 2, 3], None),
                    column("masked", vec![4, 5, 6], Some(vec![0, 1, 2])),
                ],
            }],
        }],
    };
    let expected = rmp_serde::to_vec_named(&file).expect("fixture serializes");
    let actual = consume_to_vec_named(file).expect("consuming serializer succeeds");
    assert_eq!(actual, expected);
}

fn column(name: &str, data: Vec<u8>, mask: Option<Vec<u8>>) -> EncodedColumn {
    EncodedColumn {
        name: name.to_owned(),
        data: encoded(data),
        mask: mask.map(encoded),
    }
}

fn encoded(data: Vec<u8>) -> EncodedData {
    EncodedData {
        encoding: vec![Encoding::ByteArray {
            r#type: DataType::Uint8,
        }],
        data,
    }
}
