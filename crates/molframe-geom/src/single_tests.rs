// The values compared are small integers, exact in `f32`.
#![allow(clippy::float_cmp)]

use super::*;

#[test]
fn the_products_of_the_axes_follow_the_right_hand_rule() {
    assert_eq!(cross_f32([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]), [0.0, 0.0, 1.0]);
    assert_eq!(dot_f32([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]), 32.0);
    assert_eq!(sub_f32([3.0, 2.0, 1.0], [1.0, 1.0, 1.0]), [2.0, 1.0, 0.0]);
}

#[test]
fn a_vector_too_short_or_not_finite_has_no_direction() {
    assert_eq!(normalise_f32([0.0, 3.0, 4.0], 1e-6), Some([0.0, 0.6, 0.8]));
    assert_eq!(normalise_f32([0.0, 0.0, 1e-7], 1e-6), None);
    assert_eq!(normalise_f32([f32::INFINITY, 0.0, 0.0], 1e-6), None);
}
