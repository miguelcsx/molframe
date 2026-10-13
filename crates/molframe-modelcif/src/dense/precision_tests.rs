use super::PrecisionProbe;

fn probe(values: &[f32]) -> PrecisionProbe {
    let mut probe = PrecisionProbe::new();
    for value in values {
        probe.observe(*value);
    }
    probe
}

#[test]
fn whole_numbers_need_no_decimals() {
    assert_eq!(probe(&[0.0, 7.0, 31.0]).exact_decimals(255), Some(0));
}

#[test]
fn two_decimal_scores_report_two_decimals() {
    let found = probe(&[65.38, 98.5, 91.12]);
    assert_eq!(found.exact_decimals(65_535), Some(2));
    assert!(!found.half_is_exact());
}

#[test]
fn a_ceiling_rejects_precisions_whose_integers_are_too_large() {
    assert_eq!(probe(&[31.75]).exact_decimals(255), None);
    assert_eq!(probe(&[31.75]).exact_decimals(65_535), Some(2));
}

#[test]
fn non_finite_values_are_flagged() {
    assert!(probe(&[1.0, f32::NAN]).saw_non_finite());
}

#[test]
fn dyadic_values_are_exact_in_binary16() {
    assert!(probe(&[0.5, 1.25, 31.75]).half_is_exact());
}
