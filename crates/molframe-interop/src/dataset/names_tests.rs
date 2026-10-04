use super::*;

#[test]
fn a_strategy_takes_exactly_its_own_parameter() {
    assert_eq!(
        SplitStrategy::from_parts("sequence_identity", Some(0.3), None),
        Ok(SplitStrategy::SequenceIdentity { threshold: 0.3 })
    );
    assert_eq!(
        SplitStrategy::from_parts("random", None, Some(7)),
        Ok(SplitStrategy::Random { seed: 7 })
    );
    assert_eq!(
        SplitStrategy::from_parts("temporal", None, None),
        Ok(SplitStrategy::Temporal)
    );
    assert!(
        SplitStrategy::from_parts("random", None, None).is_err(),
        "no invented seed"
    );
    assert!(SplitStrategy::from_parts("temporal", None, Some(1)).is_err());
    assert!(SplitStrategy::from_parts("sequence-identity", None, None).is_err());
    assert!(SplitStrategy::from_parts("clustered", None, None).is_err());
}
