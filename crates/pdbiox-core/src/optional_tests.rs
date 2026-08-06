use super::*;

#[test]
fn an_optional_identifier_is_no_wider_than_a_present_one() {
    assert_eq!(size_of::<OptionalSymbol>(), size_of::<SymbolId>());
    assert_eq!(size_of::<OptionalI32>(), size_of::<i32>());
}

#[test]
fn absence_and_presence_round_trip_through_the_standard_option() {
    let symbol = SymbolId::from_raw(12);
    assert_eq!(
        Option::from(OptionalSymbol::from(Some(symbol))),
        Some(symbol)
    );
    assert_eq!(Option::<SymbolId>::from(OptionalSymbol::from(None)), None);
    assert_eq!(Option::from(OptionalI32::from(Some(-3))), Some(-3));
    assert_eq!(Option::<i32>::from(OptionalI32::from(None)), None);
}

#[test]
fn negative_and_zero_sequence_positions_are_ordinary_present_values() {
    for value in [-9999, -1, 0, 1, 9999, i32::MAX] {
        assert_eq!(OptionalI32::some(value).get(), Some(value));
        assert!(OptionalI32::some(value).is_some());
    }
}

#[test]
fn the_default_is_absent_so_a_column_grown_by_resizing_reads_correctly() {
    assert!(OptionalSymbol::default().is_none());
    assert!(OptionalI32::default().is_none());
}
