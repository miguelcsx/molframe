use super::*;

#[test]
fn the_same_seed_produces_the_same_sequence() {
    let mut left = Seed::new(42);
    let mut right = Seed::new(42);
    for _ in 0..64 {
        assert_eq!(left.next_bits(), right.next_bits());
    }
}

#[test]
fn different_seeds_produce_different_sequences() {
    let mut left = Seed::new(1);
    let mut right = Seed::new(2);
    assert_ne!(left.next_bits(), right.next_bits());
}

#[test]
fn unit_values_stay_within_the_half_open_unit_interval() {
    let mut seed = Seed::new(9);
    for _ in 0..4096 {
        let value = seed.next_unit();
        assert!((0.0..1.0).contains(&value), "out of range: {value}");
    }
}

#[test]
fn signed_values_stay_within_the_requested_half_width() {
    let mut seed = Seed::new(11);
    for _ in 0..4096 {
        let value = seed.next_signed(0.25);
        assert!(value >= -0.25, "below range: {value}");
        assert!(value < 0.25, "above range: {value}");
    }
}

#[test]
fn a_derived_sequence_does_not_depend_on_how_far_the_parent_ran() {
    let parent = Seed::new(5);
    let early = parent.derive(1000).next_bits();

    let mut advanced = Seed::new(5);
    for _ in 0..100 {
        let _skipped = advanced.next_bits();
    }
    assert_eq!(Seed::new(5).derive(1000).next_bits(), early);
    assert_eq!(parent.derive(1000).next_bits(), early);
}

#[test]
fn distinct_ordinals_derive_distinct_sequences() {
    let parent = Seed::new(3);
    assert_ne!(parent.derive(0).next_bits(), parent.derive(1).next_bits());
}
