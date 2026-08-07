use super::*;

fn expression(source: &str) -> Expr {
    match parse(source) {
        Ok(parsed) => parsed.expr,
        Err(findings) => panic!("parse failed: {findings:?}"),
    }
}

#[test]
fn and_binds_more_tightly_than_or() {
    let parsed = expression("all or none and all");
    assert!(matches!(
        parsed,
        Expr::Or(_, right) if matches!(*right, Expr::And(_, _))
    ));
}

#[test]
fn a_mixed_unparenthesised_boolean_expression_warns() {
    let parsed = match parse("all or none and all") {
        Ok(parsed) => parsed,
        Err(findings) => panic!("parse failed: {findings:?}"),
    };
    assert!(
        parsed
            .warnings
            .iter()
            .any(|warning| warning.code() == Code::W4001)
    );
    let parenthesised = match parse("all or (none and all)") {
        Ok(parsed) => parsed,
        Err(findings) => panic!("parse failed: {findings:?}"),
    };
    assert!(
        !parenthesised
            .warnings
            .iter()
            .any(|warning| warning.code() == Code::W4001)
    );
}

#[test]
fn geometric_modifiers_build_typed_nodes() {
    assert!(matches!(
        expression("within 5 of (resname ATP)"),
        Expr::Geometric(GeometricExpr::Within { radius, .. }) if (radius - 5.0).abs() < f32::EPSILON
    ));
    assert!(matches!(
        expression("point 1 2 3 4"),
        Expr::Geometric(GeometricExpr::Point { point: [x, y, z], radius })
            if (x - 1.0).abs() < f32::EPSILON
                && (y - 2.0).abs() < f32::EPSILON
                && (z - 3.0).abs() < f32::EPSILON
                && (radius - 4.0).abs() < f32::EPSILON
    ));
}

#[test]
fn membership_comparison_range_and_escaped_keyword_parse() {
    assert!(matches!(
        expression("name CA CB"),
        Expr::Membership { column: Column::AtomName, values } if values.len() == 2
    ));
    assert!(matches!(
        expression("occupancy >= 0.5"),
        Expr::Comparison {
            column: Column::Occupancy,
            operator: Operator::GreaterEqual,
            ..
        }
    ));
    assert!(matches!(
        expression("resid 10:20"),
        Expr::Membership { column: Column::ResidueId, values } if values[0].as_ref() == "10:20"
    ));
    assert!(matches!(
        expression("resname \\water"),
        Expr::Membership { values, .. } if values[0].as_ref() == "water"
    ));
}

#[test]
fn flipped_comparisons_are_normalised() {
    assert!(matches!(
        expression("prop 5 > x"),
        Expr::Comparison { column: Column::X, operator: Operator::Less, value, .. }
            if (value - 5.0).abs() < f64::EPSILON
    ));
}

#[test]
fn malformed_and_unknown_selections_return_registered_diagnostics() {
    assert_eq!(
        parse("within nope of all")
            .err()
            .and_then(|findings| findings.first().map(Diagnostic::code)),
        Some(Code::E4002)
    );
    assert_eq!(
        parse("wibble A")
            .err()
            .and_then(|findings| findings.first().map(Diagnostic::code)),
        Some(Code::E4004)
    );
}
