use super::*;
use crate::{PeriodicBox, SpatialSearchOptions};
use pdbiox_core::{
    AtomSelection, ExecutionContext, MemoryBudget, ScratchPolicy, structure::UnitCell,
};

#[test]
fn counts_match_the_reference_across_workers_backends_and_periodic_images() {
    let positions: Vec<_> = (0_u16..768)
        .map(|i| [f32::from(i % 24), f32::from(i / 24), 0.0])
        .collect();
    let left = AtomSelection::range(0..600);
    let right = AtomSelection::range(100..768);
    let periodic = PeriodicBox::from_cell(UnitCell {
        lengths: [24.0, 32.0, 10.0],
        angles: [90.0; 3],
    })
    .expect("box");
    for periodic in [None, Some(&periodic)] {
        for right in [&left, &right] {
            let mut expected = None;
            for backend in [
                SpatialBackend::BruteForce,
                SpatialBackend::CellList,
                SpatialBackend::KdTree,
                SpatialBackend::NeighborList,
            ] {
                for workers in [1, 2, 4, 8] {
                    let context = ExecutionContext::builder()
                        .worker_budget(workers)
                        .build()
                        .expect("context");
                    let query = PairQuery {
                        positions: &positions,
                        left: &left,
                        right,
                        cutoff: 1.1,
                        options: SpatialSearchOptions::with_backend(backend),
                        periodic,
                        context: &context,
                    };
                    let count = count_pairs_within(&query).expect("count");
                    if let Some(expected) = expected {
                        assert_eq!(count, expected, "{backend:?}, {workers}");
                    } else {
                        expected = Some(count);
                    }
                    assert_eq!(context.reserved_bytes(), 0);
                }
            }
        }
    }
}

#[test]
fn dense_contact_output_does_not_require_pair_proportional_memory() {
    let context = ExecutionContext::builder()
        .worker_budget(8)
        .memory_budget(MemoryBudget::new(32_768).expect("budget"))
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("context");
    let positions = vec![[0.0; 3]; 1024];
    let all = AtomSelection::All(1024);
    let query = PairQuery {
        positions: &positions,
        left: &all,
        right: &all,
        cutoff: 0.0,
        options: SpatialSearchOptions::with_backend(SpatialBackend::CellList),
        periodic: None,
        context: &context,
    };
    assert_eq!(count_pairs_within(&query).expect("count"), 1024 * 1023 / 2);
    assert!(context.peak_reserved_bytes() <= 32_768);
    assert_eq!(context.reserved_bytes(), 0);
}
