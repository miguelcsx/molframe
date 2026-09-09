use super::*;
use proptest::prelude::*;
use std::sync::Arc;

#[test]
fn typed_array_codes_are_messagepack_integers_for_external_interoperability() {
    assert_eq!(
        rmp_serde::to_vec(&DataType::Int32).expect("type encodes"),
        vec![3]
    );
    assert_eq!(
        rmp_serde::from_slice::<DataType>(&[6]).expect("type decodes"),
        DataType::Uint32
    );
}

fn int32(values: &[i32]) -> Vec<u8> {
    values
        .iter()
        .flat_map(|value| value.to_le_bytes())
        .collect()
}

#[test]
fn every_numeric_inverse_in_a_mixed_chain_runs_in_reverse_order() {
    let encoded = EncodedData {
        encoding: vec![
            Encoding::FixedPoint {
                factor: 10.0,
                src_type: DataType::Float64,
            },
            Encoding::Delta {
                origin: 10,
                src_type: DataType::Int32,
            },
            Encoding::RunLength {
                src_type: DataType::Int32,
                src_size: 3,
            },
            Encoding::ByteArray {
                r#type: DataType::Int32,
            },
        ],
        data: int32(&[0, 1, 5, 2]),
    };
    assert_eq!(
        decode(&encoded).expect("chain decodes"),
        Decoded::Floats(vec![1.0, 1.5, 2.0])
    );
}

#[test]
fn string_arrays_decode_dictionary_offsets_indices_and_empty_sentinel() {
    let encoded = EncodedData {
        encoding: vec![Encoding::StringArray {
            data_encoding: vec![Encoding::ByteArray {
                r#type: DataType::Int32,
            }],
            string_data: "aAB".to_owned(),
            offset_encoding: vec![Encoding::ByteArray {
                r#type: DataType::Int32,
            }],
            offsets: int32(&[0, 1, 3]),
        }],
        data: int32(&[0, 1, -1, 0]),
    };
    let Decoded::Strings(values) = decode(&encoded).expect("strings decode") else {
        panic!("expected strings")
    };
    assert_eq!(values.iter().collect::<Vec<_>>(), ["a", "AB", "", "a"]);
    assert_eq!(values.indices(), [0, 1, 2, 0]);
    assert_eq!(values.dictionary().len(), 3);
}

#[test]
fn repeated_string_rows_use_one_shared_value_and_four_byte_indices() {
    const ROWS: usize = 100_000;

    let encoded = EncodedData {
        encoding: vec![Encoding::StringArray {
            data_encoding: vec![Encoding::ByteArray {
                r#type: DataType::Int32,
            }],
            string_data: "ALA".to_owned(),
            offset_encoding: vec![Encoding::ByteArray {
                r#type: DataType::Int32,
            }],
            offsets: int32(&[0, 3]),
        }],
        data: int32(&vec![0; ROWS]),
    };
    let Decoded::Strings(values) = decode(&encoded).expect("strings decode") else {
        panic!("expected strings")
    };

    assert_eq!(values.len(), ROWS);
    assert_eq!(values.dictionary().len(), 1);
    assert_eq!(Arc::strong_count(&values.dictionary()[0]), 1);
    assert_eq!(
        std::mem::size_of_val(values.indices()),
        ROWS * std::mem::size_of::<u32>()
    );
    assert!(values.iter().all(|value| value == "ALA"));
}

#[test]
fn public_string_column_canonicalises_dictionary_and_exposes_rows() {
    let dictionary = vec![
        Arc::<str>::from("ALA"),
        Arc::<str>::from("GLY"),
        Arc::<str>::from("ALA"),
    ];
    let values =
        DecodedStringColumn::new(dictionary, vec![2, 1, 0]).expect("valid dictionary indices");

    assert_eq!(values.dictionary().len(), 2);
    assert_eq!(values.indices(), [0, 1, 0]);
    assert_eq!(values.get(1), Some("GLY"));
    assert_eq!(
        values.iter().rev().collect::<Vec<_>>(),
        ["ALA", "GLY", "ALA"]
    );
}

#[test]
fn inconsistent_lengths_are_registered_diagnostics() {
    let encoded = EncodedData {
        encoding: vec![Encoding::ByteArray {
            r#type: DataType::Int32,
        }],
        data: vec![1, 2, 3],
    };
    assert_eq!(
        decode(&encoded).expect_err("length is bad").code(),
        Code::E1401
    );
}

proptest! {
    #[test]
    fn delta_round_trips_generated_integer_series(values in prop::collection::vec(-1000i32..1000, 0..200)) {
        let mut previous = 0i32;
        let deltas: Vec<_> = values.iter().map(|value| {
            let delta = value.wrapping_sub(previous);
            previous = *value;
            delta
        }).collect();
        let encoded = EncodedData {
            encoding: vec![
                Encoding::Delta { origin: 0, src_type: DataType::Int32 },
                Encoding::ByteArray { r#type: DataType::Int32 },
            ],
            data: int32(&deltas),
        };
        let expected = Decoded::Integers(values.into_iter().map(i64::from).collect());
        prop_assert_eq!(decode(&encoded).expect("generated chain decodes"), expected);
    }
}
