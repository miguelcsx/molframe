use super::PlddtBand;

#[test]
fn boundaries_belong_to_the_lower_band() {
    assert_eq!(PlddtBand::classify(90.0), Some(PlddtBand::Confident));
    assert_eq!(PlddtBand::classify(90.01), Some(PlddtBand::VeryHigh));
    assert_eq!(PlddtBand::classify(70.0), Some(PlddtBand::Low));
    assert_eq!(PlddtBand::classify(50.0), Some(PlddtBand::VeryLow));
    assert_eq!(PlddtBand::classify(50.5), Some(PlddtBand::Low));
    assert_eq!(PlddtBand::classify(0.0), Some(PlddtBand::VeryLow));
    assert_eq!(PlddtBand::classify(100.0), Some(PlddtBand::VeryHigh));
}

#[test]
fn non_finite_scores_have_no_band() {
    assert_eq!(PlddtBand::classify(f64::NAN), None);
    assert_eq!(PlddtBand::classify(f64::INFINITY), None);
}

#[test]
fn names_round_trip_and_bounds_are_consistent() {
    for band in PlddtBand::ALL {
        assert_eq!(PlddtBand::from_name(band.name()), Some(band));
        let above = band.lower_bound() + 0.5;
        if band != PlddtBand::VeryLow {
            assert_eq!(PlddtBand::classify(above), Some(band));
        }
    }
    assert_eq!(PlddtBand::from_name("VERY_HIGH"), Some(PlddtBand::VeryHigh));
    assert_eq!(PlddtBand::from_name("none"), None);
}

#[test]
fn bands_are_ordered_by_confidence() {
    assert!(PlddtBand::VeryLow < PlddtBand::VeryHigh);
}
