//! Corrupt, truncated and hostile BinaryCIF inputs must error, never panic or
//! allocate what they merely claim.

use molframe_bcif::{DataType, EncodedData, Encoding, decode, read, read_document};
use molframe_core::io::{InputBuffer, Limits, ReadOptions};
use proptest::prelude::*;

const SAMPLE: &[u8] = include_bytes!("../../molframe-bench/data/1crn.bcif");

fn exercise(bytes: Vec<u8>) {
    let input = InputBuffer::from_bytes(bytes);
    let _ = read(&input, &ReadOptions::default());
    if let Ok(document) = read_document(&input, Limits::default()) {
        let _ = document.to_document();
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(400))]

    #[test]
    fn truncation_never_panics(cut in 0..SAMPLE.len()) {
        exercise(SAMPLE[..cut].to_vec());
    }

    #[test]
    fn byte_flips_never_panic(
        edits in proptest::collection::vec((0..SAMPLE.len(), any::<u8>()), 1..24),
    ) {
        let mut bytes = SAMPLE.to_vec();
        for (at, value) in edits {
            bytes[at] = value;
        }
        exercise(bytes);
    }

    #[test]
    fn arbitrary_bytes_never_panic(bytes in proptest::collection::vec(any::<u8>(), 0..512)) {
        exercise(bytes);
    }
}

fn run_length_bomb(src_size: usize, count: i32) -> EncodedData {
    let mut data = Vec::new();
    data.extend_from_slice(&7_i32.to_le_bytes());
    data.extend_from_slice(&count.to_le_bytes());
    EncodedData {
        encoding: vec![
            Encoding::RunLength {
                src_type: DataType::Int32,
                src_size,
            },
            Encoding::ByteArray {
                r#type: DataType::Int32,
            },
        ],
        data,
    }
}

#[test]
fn a_run_length_declaring_an_astronomical_size_errors() {
    assert!(decode(&run_length_bomb(usize::MAX, i32::MAX)).is_err());
    assert!(decode(&run_length_bomb(1 << 40, i32::MAX)).is_err());
}

#[test]
fn a_run_length_larger_than_its_declaration_errors() {
    assert!(decode(&run_length_bomb(4, i32::MAX)).is_err());
    assert!(decode(&run_length_bomb(4, -1)).is_err());
}

#[test]
fn integer_packing_declaring_an_astronomical_size_errors() {
    let encoded = EncodedData {
        encoding: vec![
            Encoding::IntegerPacking {
                byte_count: 1,
                is_unsigned: true,
                src_size: usize::MAX,
            },
            Encoding::ByteArray {
                r#type: DataType::Uint8,
            },
        ],
        data: vec![255, 1],
    };
    assert!(decode(&encoded).is_err());
}
