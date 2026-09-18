use super::exact_integer_as_f64;

#[test]
fn binary64_boundary_is_preserved_without_rounding() {
    let boundary = 1_i64 << f64::MANTISSA_DIGITS;
    assert_eq!(
        exact_integer_as_f64(boundary),
        Some(9_007_199_254_740_992.0)
    );
    assert_eq!(
        exact_integer_as_f64(-boundary),
        Some(-9_007_199_254_740_992.0)
    );
}

#[test]
fn integers_outside_binary64_exact_domain_are_not_queried_as_neighbors() {
    let first_rounded = (1_i64 << f64::MANTISSA_DIGITS) + 1;
    assert_eq!(exact_integer_as_f64(first_rounded), None);
    assert_eq!(exact_integer_as_f64(-first_rounded), None);
}
