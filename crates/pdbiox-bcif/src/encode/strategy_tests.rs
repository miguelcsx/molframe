use super::*;
use crate::codec::{DataType, Encoding};

fn candidate(data: Vec<u8>) -> EncodedData {
    EncodedData {
        encoding: vec![Encoding::ByteArray {
            r#type: DataType::Uint8,
        }],
        data,
    }
}

#[test]
fn counting_writer_matches_exact_messagepack_length() {
    for value in [candidate(Vec::new()), candidate(vec![7; 128])] {
        let bytes = rmp_serde::to_vec(&value).expect("candidate serializes");
        assert_eq!(encoded_size(&value), bytes.len());
    }
}

#[test]
fn equal_size_candidates_keep_stable_first_choice() {
    let first = candidate(vec![1, 2, 3]);
    let second = candidate(vec![4, 5, 6]);
    let mut choice = EncodingChoice::new(first);
    choice.consider(second);
    assert_eq!(choice.finish().data, vec![1, 2, 3]);
}

#[test]
fn strictly_smaller_candidate_replaces_the_current_choice() {
    let current = candidate(vec![1; 128]);
    let smaller = candidate(vec![2]);
    let mut choice = EncodingChoice::new(current);
    choice.consider(smaller);
    assert_eq!(choice.finish().data, vec![2]);
}

#[test]
fn payload_lower_bound_matches_exact_messagepack_size() {
    for length in [0, 255, 256, 65_535, 65_536] {
        let value = candidate(vec![7; length]);
        assert_eq!(
            encoded_size_lower_bound(value.encoding.clone(), length),
            Some(encoded_size(&value))
        );
    }
}
