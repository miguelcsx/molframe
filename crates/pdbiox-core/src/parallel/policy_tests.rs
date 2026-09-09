use super::*;

#[test]
fn the_default_policy_is_deterministic() {
    assert_eq!(ReductionPolicy::default(), ReductionPolicy::Deterministic);
    assert!(ReductionPolicy::default().is_deterministic());
}

#[test]
fn the_fast_policy_is_not_deterministic() {
    assert!(!ReductionPolicy::Fast.is_deterministic());
}

#[test]
fn each_policy_has_a_stable_recorded_name() {
    assert_eq!(ReductionPolicy::Deterministic.name(), "deterministic");
    assert_eq!(ReductionPolicy::Fast.name(), "fast");
    assert_eq!(ReductionPolicy::Fast.to_string(), "fast");
}
