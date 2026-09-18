use super::*;

#[test]
fn ranges_lists_and_juxtaposition_expand_as_an_ordered_cartesian_product() {
    let expression = match OperExpression::parse("(1-3,5)(A,B)") {
        Ok(expression) => expression,
        Err(error) => panic!("expression failed: {error}"),
    };
    let mut combinations = Vec::new();
    expression.for_each_combination(|values| combinations.push(values.join("/")));
    assert_eq!(
        combinations,
        ["1/A", "1/B", "2/A", "2/B", "3/A", "3/B", "5/A", "5/B"]
    );
    assert_eq!(expression.combination_count(), 8);
}

#[test]
fn one_unparenthesised_list_is_one_factor() {
    let expression = OperExpression::parse("1,3-4").expect("valid expression");
    assert_eq!(expression.combination_count(), 3);
    assert_eq!(expression.factors()[0].len(), 3);
}

#[test]
fn malformed_or_oversized_products_are_diagnostics_before_expansion() {
    assert_eq!(
        OperExpression::parse("(4-1)")
            .err()
            .map(|error| error.code()),
        Some(Code::E6010)
    );
    assert_eq!(
        OperExpression::parse_with_limit("(1-60)(61-88)", 1_000)
            .err()
            .map(|error| error.code()),
        Some(Code::E6011)
    );
}

#[test]
fn dictionary_character_codes_are_not_restricted_to_alphanumerics() {
    let expression = match OperExpression::parse("(X0-foo,a.b:1)") {
        Ok(expression) => expression,
        Err(error) => panic!("character-code expression failed: {error}"),
    };
    let values = expression.factors()[0]
        .iter()
        .map(Box::as_ref)
        .collect::<Vec<_>>();
    assert_eq!(values, ["X0-foo", "a.b:1"]);
}
