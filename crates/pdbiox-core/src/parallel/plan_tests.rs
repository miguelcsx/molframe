use super::*;

#[test]
fn blocks_cover_every_item_exactly_once() {
    for count in [0_usize, 1, 7, 64, 65, 1000] {
        for block in [1_usize, 3, 64, 128] {
            let plan = BlockPlan::new(count, block);
            let mut seen = 0;
            for range in plan.ranges() {
                assert_eq!(range.start, seen, "gap or overlap at {seen}");
                seen = range.end;
            }
            assert_eq!(seen, count, "count {count} block {block}");
        }
    }
}

#[test]
fn the_partition_does_not_depend_on_a_worker_count() {
    let plan = BlockPlan::new(1000, 64);
    let ranges: Vec<_> = plan.ranges().collect();
    for workers in [1_usize, 2, 4, 16, 64] {
        let _useful = plan.useful_workers(workers);
        assert_eq!(plan.ranges().collect::<Vec<_>>(), ranges);
    }
}

#[test]
fn a_zero_block_size_is_raised_to_one() {
    let plan = BlockPlan::new(3, 0);
    assert_eq!(plan.block(), 1);
    assert_eq!(plan.blocks(), 3);
}

#[test]
fn an_empty_plan_has_no_blocks_and_no_ranges() {
    let plan = BlockPlan::new(0, 8);
    assert!(plan.is_empty());
    assert_eq!(plan.blocks(), 0);
    assert_eq!(plan.ranges().count(), 0);
    assert_eq!(plan.range(0), None);
}

#[test]
fn a_block_past_the_end_has_no_range() {
    let plan = BlockPlan::new(10, 4);
    assert_eq!(plan.range(2), Some(8..10));
    assert_eq!(plan.range(3), None);
}

#[test]
fn workers_never_exceed_the_block_count() {
    let plan = BlockPlan::new(10, 4);
    assert_eq!(plan.useful_workers(64), 3);
    assert_eq!(plan.useful_workers(0), 1);
    assert_eq!(plan.useful_workers(2), 2);
}

#[test]
fn the_default_block_is_the_documented_leaf_size() {
    assert_eq!(
        BlockPlan::with_default_block(128).block(),
        DEFAULT_BLOCK_ITEMS
    );
}
