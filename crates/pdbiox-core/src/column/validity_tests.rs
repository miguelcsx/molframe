use super::*;

#[test]
fn an_all_present_mask_costs_a_length_rather_than_two_bit_sets() {
    let mask = ValidityMask::all_present(1_000_000);
    assert!(mask.is_all_present());
    assert_eq!(mask.present_count(), 1_000_000);
    assert!(matches!(mask, ValidityMask::AllPresent(_)));
}

#[test]
fn unrecorded_and_inapplicable_are_distinguishable_after_a_round_trip() {
    let mask: ValidityMask = [Presence::Present, Presence::Unknown, Presence::Inapplicable]
        .into_iter()
        .collect();
    assert_eq!(mask.get(0), Presence::Present);
    assert_eq!(mask.get(1), Presence::Unknown);
    assert_eq!(mask.get(2), Presence::Inapplicable);
    assert_eq!(mask.present_count(), 1);
}

#[test]
fn marking_a_position_absent_materialises_the_bit_sets_once() {
    let mut mask = ValidityMask::all_present(4);
    mask.set(2, Presence::Unknown);
    assert!(matches!(mask, ValidityMask::Mixed { .. }));
    assert_eq!(mask.get(0), Presence::Present);
    assert_eq!(mask.get(2), Presence::Unknown);
    assert_eq!(mask.present_count(), 3);
}

#[test]
fn marking_everything_present_again_collapses_the_mask() {
    let mut mask = ValidityMask::all_present(4);
    mask.set(1, Presence::Inapplicable);
    mask.set(1, Presence::Present);
    mask.compact();
    assert!(matches!(mask, ValidityMask::AllPresent(4)));
}

#[test]
fn a_position_past_the_end_is_inapplicable_rather_than_present() {
    let mask = ValidityMask::all_present(2);
    assert_eq!(mask.get(2), Presence::Inapplicable);
    assert_eq!(mask.get(u32::MAX), Presence::Inapplicable);
}

#[test]
fn setting_present_on_a_compact_mask_leaves_it_compact() {
    let mut mask = ValidityMask::all_present(8);
    mask.set(3, Presence::Present);
    assert!(matches!(mask, ValidityMask::AllPresent(8)));
}
