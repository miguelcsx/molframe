use super::*;
use crate::codec::{Decoded, decode};
use proptest::prelude::*;

#[test]
fn strings_round_trip_in_first_seen_dictionary_order() {
    let values = vec!["ALA".to_owned(), "GLY".to_owned(), "ALA".to_owned()];
    let encoded = encode_strings(&values).expect("values encode");
    assert_eq!(
        decode(&encoded).expect("values decode"),
        Decoded::Strings(values)
    );
}

proptest! {
    #[test]
    fn every_string_column_round_trips(values in proptest::collection::vec("[A-Za-z0-9_]{0,12}", 0..128)) {
        let encoded = encode_strings(&values).expect("values encode");
        prop_assert_eq!(decode(&encoded).expect("values decode"), Decoded::Strings(values));
    }
}
