#![allow(
    clippy::float_cmp,
    reason = "these read back stored values, not computed ones"
)]

use super::*;

#[test]
fn a_lane_is_three_cache_lines_so_the_allocation_is_aligned() {
    assert_eq!(size_of::<CoordLane>(), 192);
    assert_eq!(align_of::<CoordLane>(), 64);
}

#[test]
fn positions_are_contiguous_across_lane_boundaries() {
    let block: CoordinateBlock = (0..40).map(|i| [i as f32, 0.0, 0.0]).collect();
    assert_eq!(block.len(), 40);
    let slice = block.as_slice();
    assert_eq!(slice.len(), 40);
    for (i, position) in slice.iter().enumerate() {
        assert_eq!(position[0], i as f32);
    }
}

#[test]
fn the_slice_stops_at_the_position_count_not_at_the_lane_boundary() {
    let mut block = CoordinateBlock::new();
    block.push([1.0, 1.0, 1.0]);
    assert_eq!(block.as_slice().len(), 1);
    assert_eq!(block.allocated_bytes(), 192);
}

#[test]
fn a_cloned_block_shares_storage_until_one_copy_is_edited() {
    let original: CoordinateBlock = [[1.0, 2.0, 3.0]].into_iter().collect();
    let mut edited = original.clone();
    assert!(std::ptr::eq(
        original.as_slice().as_ptr(),
        edited.as_slice().as_ptr()
    ));

    edited.as_mut_slice()[0][0] = 9.0;
    assert!(!std::ptr::eq(
        original.as_slice().as_ptr(),
        edited.as_slice().as_ptr()
    ));
    assert_eq!(original.as_slice()[0][0], 1.0);
    assert_eq!(edited.as_slice()[0][0], 9.0);
}

#[test]
fn a_range_past_the_end_yields_nothing_rather_than_a_short_slice() {
    let block: CoordinateBlock = (0..5).map(|_| [0.0; 3]).collect();
    assert!(block.range(0..5).is_some());
    assert!(block.range(0..6).is_none());
    assert!(block.range(3..99).is_none());
}

#[test]
fn a_bounding_box_ignores_non_finite_positions() {
    let mut bounds = Aabb::EMPTY;
    bounds.extend([1.0, 2.0, 3.0]);
    bounds.extend([f32::NAN, -1.0, 0.0]);
    assert_eq!(bounds.min, [1.0, -1.0, 0.0]);
    assert_eq!(bounds.max, [1.0, 2.0, 3.0]);
}

#[test]
fn an_empty_bounding_box_is_never_within_reach_of_anything() {
    let full = Aabb {
        min: [0.0; 3],
        max: [1.0; 3],
    };
    assert!(!Aabb::EMPTY.within(&full, 1000.0));
    assert!(!full.within(&Aabb::EMPTY, 1000.0));
}

#[test]
fn boxes_that_are_further_apart_than_the_cutoff_are_rejected() {
    let left = Aabb {
        min: [0.0; 3],
        max: [1.0; 3],
    };
    let right = Aabb {
        min: [10.0, 0.0, 0.0],
        max: [11.0, 1.0, 1.0],
    };
    assert!(!left.within(&right, 8.0));
    assert!(left.within(&right, 9.0));
    assert!(left.within(&left, 0.0));
}

#[test]
fn a_generation_advances_and_stalls_rather_than_wrapping() {
    assert_eq!(CoordinateGeneration::INITIAL.next().get(), 1);
    assert_eq!(
        CoordinateGeneration::INITIAL.next().next(),
        CoordinateGeneration::INITIAL.next().next()
    );
}
