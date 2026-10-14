//! Pairwise bonding predicates and per-element reach used by bond perception.
//!
//! Every test here is constant time and works on squared distances, so the
//! grid search stays free of square roots; the thresholds are reference data
//! whose per-element reach is a [`crate::BondReach`] table.

use molframe_core::structure::AtomRef;

/// Whether a distance lies between the closest plausible contact and the sum of
/// the two elements' bonding thresholds.
pub(super) fn in_bonding_window(
    left_threshold: f32,
    right_threshold: f32,
    distance_squared: f32,
) -> bool {
    let maximum = maximum_bonding_distance(left_threshold, right_threshold);
    distance_squared >= 0.16 && distance_squared <= maximum * maximum
}

pub(super) fn maximum_bonding_distance(left_threshold: f32, right_threshold: f32) -> f32 {
    (left_threshold + right_threshold) / 1.95
}

/// Whether two atoms may bond at all: two hydrogens never do, and two atoms in
/// different alternate locations do not coexist.
pub(super) fn pair_compatible(left: AtomRef<'_>, right: AtomRef<'_>) -> bool {
    let both_hydrogen = left
        .element()
        .is_some_and(molframe_core::Element::is_hydrogen)
        && right
            .element()
            .is_some_and(molframe_core::Element::is_hydrogen);
    !both_hydrogen && alt_compatible(left, right)
}

pub(super) fn alt_compatible(left: AtomRef<'_>, right: AtomRef<'_>) -> bool {
    match (left.alt_id(), right.alt_id()) {
        (Some(left), Some(right)) => left.is_blank() || right.is_blank() || left == right,
        _ => true,
    }
}

pub(super) fn within(left: AtomRef<'_>, right: AtomRef<'_>, maximum: f32) -> bool {
    let (Some(left), Some(right)) = (left.position(), right.position()) else {
        return false;
    };
    squared_distance(left, right) <= maximum * maximum
}

pub(super) fn squared_distance(left: [f32; 3], right: [f32; 3]) -> f32 {
    left.iter().zip(right).map(|(a, b)| (a - b) * (a - b)).sum()
}
