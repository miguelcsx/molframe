use super::*;
use crate::codec::{Decoded, decode};
use proptest::prelude::*;

#[test]
fn strings_round_trip_in_first_seen_dictionary_order() {
    let values = vec!["ALA".to_owned(), "GLY".to_owned(), "ALA".to_owned()];
    let encoded = encode_strings(&values).expect("values encode");
    let Decoded::Strings(decoded) = decode(&encoded).expect("values decode") else {
        panic!("expected strings")
    };
    assert!(decoded.iter().eq(values.iter().map(String::as_str)));
    assert_eq!(decoded.dictionary().len(), 2);
}

proptest! {
    #[test]
    fn every_string_column_round_trips(values in proptest::collection::vec("[A-Za-z0-9_]{0,12}", 0..128)) {
        let encoded = encode_strings(&values).expect("values encode");
        let Decoded::Strings(decoded) = decode(&encoded).expect("values decode") else {
            return Err(TestCaseError::fail("expected strings"));
        };
        prop_assert!(decoded.iter().eq(values.iter().map(String::as_str)));
    }
}
