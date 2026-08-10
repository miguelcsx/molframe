use crate::predicate_pattern::{NumericPattern, parse_residue, split_range};

#[test]
fn signed_numbers_and_insertion_codes_split_without_losing_the_sign() {
    let negative = match parse_residue("-12A") {
        Ok(value) => value,
        Err(error) => panic!("parse failed: {error}"),
    };
    assert_eq!(negative.number, -12);
    assert_eq!(negative.insertion.as_ref(), "A");
    assert_eq!(split_range("-12A--10B"), Some(("-12A", "-10B")));
}

#[test]
fn numeric_ranges_are_inclusive_and_order_independent() {
    let pattern = match NumericPattern::parse("5:2") {
        Ok(pattern) => pattern,
        Err(error) => panic!("parse failed: {error}"),
    };
    assert!(pattern.matches(2.0));
    assert!(pattern.matches(5.0));
    assert!(!pattern.matches(5.1));
}
