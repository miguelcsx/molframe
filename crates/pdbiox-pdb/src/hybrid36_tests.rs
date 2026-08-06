use super::*;

#[test]
fn ordinary_decimal_values_read_and_write_unchanged() {
    assert_eq!(decode("    1", 5), Some(1));
    assert_eq!(decode("99999", 5), Some(99_999));
    assert_eq!(encode(1, 5).as_deref(), Some("    1"));
    assert_eq!(encode(99_999, 5).as_deref(), Some("99999"));
}

#[test]
fn the_sequence_continues_into_upper_case_where_the_decimal_range_stops() {
    assert_eq!(encode(100_000, 5).as_deref(), Some("A0000"));
    assert_eq!(decode("A0000", 5), Some(100_000));
    assert_eq!(encode(100_001, 5).as_deref(), Some("A0001"));
    assert_eq!(decode("A0001", 5), Some(100_001));
}

#[test]
fn residue_numbers_use_the_same_scheme_in_four_columns() {
    assert_eq!(encode(9_999, 4).as_deref(), Some("9999"));
    assert_eq!(encode(10_000, 4).as_deref(), Some("A000"));
    assert_eq!(decode("A000", 4), Some(10_000));
    assert_eq!(decode("9999", 4), Some(9_999));
}

#[test]
fn every_value_across_both_case_ranges_survives_a_round_trip() {
    let upper_span = 36i64.pow(5) - 10 * 36i64.pow(4);
    for value in [
        100_000,
        100_001,
        123_456,
        1_000_000,
        100_000 + upper_span - 1,
    ] {
        let encoded = encode(value, 5);
        let decoded = encoded.as_deref().and_then(|text| decode(text, 5));
        assert_eq!(decoded, Some(value), "{value} encoded as {encoded:?}");
    }
    // The first value of the lower-case range, and one beyond it.
    for value in [100_000 + upper_span, 100_000 + upper_span + 1] {
        let encoded = encode(value, 5);
        assert!(
            encoded
                .as_deref()
                .is_some_and(|text| text.bytes().any(|b| b.is_ascii_lowercase())),
            "{value} should have reached the lower-case range, got {encoded:?}"
        );
        assert_eq!(
            encoded.as_deref().and_then(|text| decode(text, 5)),
            Some(value)
        );
    }
}

#[test]
fn negative_numbers_stay_decimal_because_the_scheme_has_no_place_for_them() {
    assert_eq!(encode(-12, 5).as_deref(), Some("  -12"));
    assert_eq!(decode("  -12", 5), Some(-12));
    assert_eq!(decode("-999", 4), Some(-999));
}

#[test]
fn a_field_that_cannot_be_read_yields_nothing_rather_than_a_wrong_number() {
    assert_eq!(decode("", 5), None);
    assert_eq!(decode("   ", 5), None);
    assert_eq!(decode("1!2", 5), None);
}

#[test]
fn the_encoder_says_when_a_value_still_will_not_fit() {
    assert!(encode(i64::MAX, 5).is_none());
    assert!(encode(-1_000_000, 5).is_none());
}

#[test]
fn a_value_inside_the_decimal_range_never_needs_the_extended_scheme() {
    assert!(!needs_hybrid36(99_999, 5));
    assert!(needs_hybrid36(100_000, 5));
    assert!(!needs_hybrid36(9_999, 4));
    assert!(needs_hybrid36(10_000, 4));
}
