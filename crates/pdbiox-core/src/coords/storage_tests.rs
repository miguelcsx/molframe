use super::*;
use num_traits::ToPrimitive;

#[test]
fn a_lane_is_three_cache_lines_so_the_allocation_is_aligned() {
    assert_eq!(size_of::<CoordLane>(), 192);
    assert_eq!(align_of::<CoordLane>(), 64);
}

#[test]
fn positions_are_contiguous_across_lane_boundaries() {
    let block: CoordinateBlock = (0..40)
        .map(|i| [i.to_f32().expect("small coordinate"), 0.0, 0.0])
        .collect();
    assert_eq!(block.len(), 40);
    let slice = block.as_slice();
    assert_eq!(slice.len(), 40);
    for (i, position) in slice.iter().enumerate() {
        assert_eq!(
            position[0].to_bits(),
            i.to_f32().expect("small coordinate").to_bits()
        );
    }
}

#[test]
fn the_slice_stops_at_the_position_count_not_at_the_lane_boundary() {
    let mut block = CoordinateBlock::new();
    block.push([1.0, 1.0, 1.0]);
    assert_eq!(block.as_slice().len(), 1);
    assert_eq!(block.allocated_bytes(), block.lanes.capacity() * 192);
}

#[test]
fn allocated_bytes_include_reserved_lanes_before_positions_are_written() {
    let mut block = CoordinateBlock::with_capacity(1000);
    let bytes = block.allocated_bytes();
    assert!(bytes >= 1000 * size_of::<[f32; 3]>());
    assert!(block.is_empty());
    block.push([1.0; 3]);
    assert_eq!(block.allocated_bytes(), bytes);
    let shared = block.clone();
    assert_eq!(shared.allocated_bytes(), bytes);
    drop(block);
    assert_eq!(shared.allocated_bytes(), bytes);
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
    assert_eq!(original.as_slice()[0][0].to_bits(), 1.0_f32.to_bits());
    assert_eq!(edited.as_slice()[0][0].to_bits(), 9.0_f32.to_bits());
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
    assert_eq!(
        bounds.min.map(f32::to_bits),
        [1.0, -1.0, 0.0].map(f32::to_bits)
    );
    assert_eq!(
        bounds.max.map(f32::to_bits),
        [1.0, 2.0, 3.0].map(f32::to_bits)
    );
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
fn a_generation_advances_and_rejects_overflow() {
    assert_eq!(
        CoordinateGeneration::INITIAL
            .next()
            .map(CoordinateGeneration::get),
        Some(1)
    );
    assert_eq!(
        CoordinateGeneration::INITIAL
            .next()
            .and_then(CoordinateGeneration::next),
        Some(CoordinateGeneration(2))
    );
    assert_eq!(CoordinateGeneration(u64::MAX).next(), None);
}
