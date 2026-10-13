#![allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]

use super::{crowding, knn_within};
use crate::{KnnOptions, KnnResult, PeriodicBox, SpatialBackend, pairs_within};
use molframe_core::structure::UnitCell;
use molframe_core::{AtomSelection, ExecutionContext};

fn points(count: usize, seed: u64, scale: f32) -> Vec<[f32; 3]> {
    let mut state = seed;
    let mut next = move || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((state >> 40) as f32) / ((1_u64 << 24) as f32) * scale
    };
    (0..count).map(|_| [next(), next(), next()]).collect()
}

fn run(
    positions: &[[f32; 3]],
    query: &AtomSelection,
    targets: &AtomSelection,
    options: KnnOptions,
    backend: SpatialBackend,
    periodic: Option<&PeriodicBox>,
) -> KnnResult {
    knn_within(
        positions,
        query,
        targets,
        5,
        options,
        backend,
        periodic,
        &ExecutionContext::default(),
    )
    .expect("knn")
}

#[test]
fn tree_and_brute_force_agree_for_selections_and_limits() {
    let positions = points(900, 7, 20.0);
    let query = AtomSelection::from_sorted((0..900).step_by(3).collect());
    let targets = AtomSelection::from_sorted((0..900).filter(|i| i % 4 != 0).collect());
    for options in [
        KnnOptions::default(),
        KnnOptions {
            max_distance: Some(2.0),
            exclude_self: true,
        },
    ] {
        let tree = run(
            &positions,
            &query,
            &targets,
            options,
            SpatialBackend::Auto,
            None,
        );
        let brute = run(
            &positions,
            &query,
            &targets,
            options,
            SpatialBackend::BruteForce,
            None,
        );
        assert_eq!(tree.indices, brute.indices);
        assert_eq!(tree.distances_squared, brute.distances_squared);
        assert_eq!(tree.rows(), 300);
    }
}

#[test]
fn periodic_queries_use_minimum_image_distances() {
    let positions = points(400, 11, 15.0);
    let cell = PeriodicBox::from_cell(UnitCell {
        lengths: [15.0, 16.0, 17.0],
        angles: [90.0, 90.0, 90.0],
    })
    .expect("cell");
    let all = AtomSelection::All(400);
    let options = KnnOptions {
        max_distance: Some(4.0),
        exclude_self: true,
    };
    let tree = run(
        &positions,
        &all,
        &all,
        options,
        SpatialBackend::KdTree,
        Some(&cell),
    );
    let brute = run(
        &positions,
        &all,
        &all,
        options,
        SpatialBackend::BruteForce,
        Some(&cell),
    );
    assert_eq!(tree.indices, brute.indices);
    assert_eq!(tree.distances_squared, brute.distances_squared);
}

#[test]
fn invalid_distance_and_out_of_range_queries_fail() {
    let positions = points(10, 3, 5.0);
    let all = AtomSelection::All(10);
    for backend in [SpatialBackend::Auto, SpatialBackend::BruteForce] {
        let options = KnnOptions {
            max_distance: Some(f32::NAN),
            exclude_self: false,
        };
        let result = knn_within(
            &positions,
            &all,
            &all,
            3,
            options,
            backend,
            None,
            &ExecutionContext::default(),
        );
        assert!(result.is_err());
    }
}

#[test]
fn crowding_probe_flags_dense_cores_only() {
    let uniform = points(20_000, 5, 200.0);
    let atoms: Vec<u32> = (0..20_000).collect();
    assert!(!crowding::is_crowded(&uniform, &atoms, 5.0));
    let mut dense = points(20_000, 5, 200.0);
    for (index, point) in dense.iter_mut().enumerate().take(10_000) {
        *point = [index as f32 * 1.0e-4, 0.0, 0.0];
    }
    assert!(crowding::is_crowded(&dense, &atoms, 5.0));
    assert!(!crowding::is_crowded(&dense, &atoms[..100], 5.0));
}

#[test]
fn auto_on_crowded_data_matches_the_cell_list() {
    let mut positions = points(6_000, 9, 120.0);
    for (index, point) in positions.iter_mut().enumerate().take(3_000) {
        *point = [index as f32 * 2.0e-3, 0.0, 1.0];
    }
    let all = AtomSelection::All(6_000);
    let ctx = ExecutionContext::default();
    let auto = pairs_within(
        &positions,
        &all,
        &all,
        1.5,
        SpatialBackend::Auto,
        None,
        &ctx,
    )
    .expect("auto");
    let cell = pairs_within(
        &positions,
        &all,
        &all,
        1.5,
        SpatialBackend::CellList,
        None,
        &ctx,
    )
    .expect("cell");
    assert_eq!(auto, cell);
}
