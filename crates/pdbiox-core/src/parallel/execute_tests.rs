use super::*;

#[test]
fn results_arrive_in_block_order_at_every_worker_count() {
    let plan = BlockPlan::new(1000, 64);
    let serial = match map_blocks(plan, 1, |index, range| (index, range.start, range.end)) {
        Ok(serial) => serial,
        Err(error) => panic!("serial run failed: {error}"),
    };

    for workers in [2_usize, 3, 4, 8, 16, 64] {
        match map_blocks(plan, workers, |index, range| {
            (index, range.start, range.end)
        }) {
            Ok(parallel) => assert_eq!(parallel, serial, "differs at {workers} workers"),
            Err(error) => panic!("run at {workers} workers failed: {error}"),
        }
    }
}

#[test]
fn a_floating_point_reduction_is_bit_identical_at_every_worker_count() {
    // Deliberately mixed magnitudes: a reassociated sum loses the small terms,
    // so this fails if the block boundaries ever move with the worker count.
    let values: Vec<f64> = (0..4096)
        .map(|index| if index % 2 == 0 { 1.0e10 } else { 0.1 })
        .collect();
    let plan = BlockPlan::new(values.len(), 64);

    let sum_at = |workers: usize| -> f64 {
        let Ok(partials) = map_blocks(plan, workers, |_, range| {
            values[range]
                .iter()
                .fold(0.0_f64, |total, value| total + value)
        }) else {
            panic!("run at {workers} workers failed")
        };
        partials.iter().fold(0.0_f64, |total, value| total + value)
    };

    let serial = sum_at(1);
    for workers in [2_usize, 4, 8, 16] {
        assert_eq!(
            sum_at(workers).to_bits(),
            serial.to_bits(),
            "sum differs at {workers} workers"
        );
    }
}

#[test]
fn an_empty_plan_produces_nothing() {
    let plan = BlockPlan::new(0, 8);
    match map_blocks(plan, 4, |index, _| index) {
        Ok(produced) => assert!(produced.is_empty()),
        Err(error) => panic!("empty run failed: {error}"),
    }
}

#[test]
fn one_worker_visits_every_block() {
    let plan = BlockPlan::new(10, 4);
    match map_blocks(plan, 1, |index, _| index) {
        Ok(produced) => assert_eq!(produced, vec![0, 1, 2]),
        Err(error) => panic!("serial run failed: {error}"),
    }
}

#[test]
fn more_workers_than_blocks_still_visits_every_block_once() {
    let plan = BlockPlan::new(10, 4);
    match map_blocks(plan, 64, |index, _| index) {
        Ok(produced) => assert_eq!(produced, vec![0, 1, 2]),
        Err(error) => panic!("oversubscribed run failed: {error}"),
    }
}

#[test]
fn context_execution_is_deterministic_for_nested_calls() {
    let context = ExecutionContext::builder()
        .worker_budget(4)
        .build()
        .expect("execution context");
    let outer = BlockPlan::new(256, 64);
    let actual = map_blocks_in(outer, &context, |outer_index, _| {
        let inner = map_blocks_in(BlockPlan::new(16, 4), &context, |inner_index, _| {
            outer_index * 10 + inner_index
        });
        inner.expect("nested shared-pool execution")
    })
    .expect("outer shared-pool execution");
    assert_eq!(actual[0], [0, 1, 2, 3]);
    assert_eq!(actual[3], [30, 31, 32, 33]);
}
