use super::*;

#[test]
fn contiguous_cells_find_only_pairs_inside_the_cutoff() {
    let positions = [[0.0, 0.0, 0.0], [0.9, 0.0, 0.0], [4.0, 0.0, 0.0]];
    let list = match CellList::build(&positions, &[0, 1, 2], 1.0, None) {
        Ok(list) => list,
        Err(error) => panic!("build failed: {error}"),
    };
    let pairs = match list.pairs(&[0, 1, 2], 1.0) {
        Ok(pairs) => pairs,
        Err(error) => panic!("query failed: {error}"),
    };
    assert_eq!(pairs.len(), 1);
    assert_eq!((pairs[0].first, pairs[0].second), (0, 1));
}

#[test]
fn a_query_cannot_exceed_the_cutoff_the_grid_was_built_for() {
    let positions = [[0.0, 0.0, 0.0]];
    let list = match CellList::build(&positions, &[0], 1.0, None) {
        Ok(list) => list,
        Err(error) => panic!("build failed: {error}"),
    };
    assert_eq!(list.pairs(&[0], 2.0), Err(SpatialError::InvalidCutoff));
}

#[test]
fn cell_grid_coarsening_obeys_the_public_memory_policy() {
    let options = CellGridOptions {
        maximum_cell_count: 1,
        edge_growth_factor: 4.0,
    };
    let Ok((_, dimensions, count)) =
        geometry::grid_geometry([0.0, 0.0, 0.0], [3.0, 3.0, 3.0], 1.0, options)
    else {
        panic!("finite test geometry must be representable");
    };
    assert_eq!(dimensions, [1, 1, 1]);
    assert_eq!(count, 1);
}

#[test]
fn query_coordinates_below_target_bounds_clamp_to_the_first_cell() {
    let positions = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let list = match CellList::build(&positions, &[2], 1.0, None) {
        Ok(list) => list,
        Err(error) => panic!("build failed: {error}"),
    };
    let pairs = match list.pairs(&[0, 1], 1.0) {
        Ok(pairs) => pairs,
        Err(error) => panic!("query failed: {error}"),
    };
    assert_eq!(pairs.len(), 2);
}

#[test]
fn arbitrary_points_include_coincident_targets_and_obey_the_radius() {
    let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.01, 0.0, 0.0]];
    let list = CellList::build(&positions, &[0, 1, 2], 1.0, None).unwrap();
    let mut neighbors = Vec::new();
    list.for_each_neighbor([0.0; 3], 1.0, |atom, squared| {
        neighbors.push((atom, squared));
    })
    .unwrap();
    neighbors.sort_by_key(|(atom, _)| *atom);
    assert_eq!(neighbors, [(0, 0.0), (1, 1.0)]);
    neighbors.clear();
    list.for_each_neighbor([-0.5, 0.0, 0.0], 1.0, |atom, squared| {
        neighbors.push((atom, squared));
    })
    .unwrap();
    assert_eq!(neighbors, [(0, 0.25)]);
    assert_eq!(
        list.for_each_neighbor([f32::NAN; 3], 1.0, |_, _| {}),
        Err(SpatialError::NumericRangeExceeded)
    );
}

#[test]
fn arbitrary_periodic_points_include_wrapped_targets_once() {
    let positions = [[0.1, 0.0, 0.0], [5.0, 0.0, 0.0]];
    let periodic = PeriodicBox::from_cell(molframe_core::structure::UnitCell {
        lengths: [10.0; 3],
        angles: [90.0; 3],
    })
    .unwrap();
    let list = CellList::build(&positions, &[0, 1], 1.0, Some(periodic)).unwrap();
    let mut neighbors = Vec::new();
    list.for_each_neighbor([9.9, 0.0, 0.0], 1.0, |atom, squared| {
        neighbors.push((atom, squared));
    })
    .unwrap();
    assert_eq!(neighbors.len(), 1);
    assert_eq!(neighbors[0].0, 0);
    assert!((neighbors[0].1 - 0.04).abs() < 1.0e-5);
}

/// Deterministic scatter in a cube of `extent` ångström, no external RNG.
#[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]
fn scatter(count: usize, extent: f32, seed: u64) -> Vec<[f32; 3]> {
    let mut state = seed;
    let mut next = move || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((state >> 40) as f32 / (1u64 << 24) as f32) * extent
    };
    (0..count).map(|_| [next(), next(), next()]).collect()
}

#[allow(clippy::cast_possible_truncation)]
fn brute_pairs(positions: &[[f32; 3]], cutoff: f32) -> Vec<(u32, u32)> {
    let mut found = Vec::new();
    for left in 0..positions.len() {
        for right in left + 1..positions.len() {
            let squared: f32 = (0..3)
                .map(|axis| (positions[left][axis] - positions[right][axis]).powi(2))
                .sum();
            if squared <= cutoff * cutoff {
                found.push((left as u32, right as u32));
            }
        }
    }
    found
}

#[allow(clippy::cast_possible_truncation)]
fn listed_pairs(list: &CellList<'_>, count: usize, cutoff: f32) -> Vec<(u32, u32)> {
    let all: Vec<u32> = (0..count as u32).collect();
    let mut pairs: Vec<(u32, u32)> = list
        .pairs(&all, cutoff)
        .unwrap()
        .into_iter()
        .map(|pair| (pair.first.min(pair.second), pair.first.max(pair.second)))
        .collect();
    pairs.sort_unstable();
    pairs.dedup();
    pairs
}

#[test]
#[allow(clippy::cast_possible_truncation)]
fn a_far_outlier_selects_the_sparse_layout_and_changes_no_pair() {
    let mut positions = scatter(1500, 30.0, 7);
    let near = brute_pairs(&positions, 3.5);
    for outlier in [1.0e4_f32, 1.0e6, 1.0e9] {
        positions.push([outlier, -outlier, outlier]);
        positions.push([outlier + 1.0, -outlier, outlier]);
        let targets: Vec<u32> = (0..positions.len() as u32).collect();
        let list = CellList::build(&positions, &targets, 3.5, None).unwrap();
        let Some(CellGrid::Cartesian(grid)) = &list.grid else {
            panic!("a Cartesian grid");
        };
        assert!(grid.sparse, "outlier {outlier} must not force a dense grid");
        assert!(grid.offsets.len() <= positions.len() + 1);
        assert_eq!(
            listed_pairs(&list, positions.len(), 3.5),
            brute_pairs(&positions, 3.5)
        );
        positions.truncate(1500);
        assert_eq!(brute_pairs(&positions, 3.5), near);
    }
}

#[test]
#[allow(clippy::cast_possible_truncation)]
fn sparse_and_dense_layouts_agree_with_brute_force() {
    let positions = scatter(1200, 60.0, 11);
    let targets: Vec<u32> = (0..positions.len() as u32).collect();
    let dense = CellList::build(&positions, &targets, 4.0, None).unwrap();
    let Some(CellGrid::Cartesian(dense_grid)) = &dense.grid else {
        panic!("a Cartesian grid");
    };
    assert!(!dense_grid.sparse);
    // Zero dense allowance forces the sparse layout over the same atoms.
    let forced = CellGridOptions {
        maximum_cell_count: 1,
        edge_growth_factor: 2.0,
    };
    let sparse = CellList::build_with_options(&positions, &targets, 4.0, None, forced).unwrap();
    let Some(CellGrid::Cartesian(sparse_grid)) = &sparse.grid else {
        panic!("a Cartesian grid");
    };
    assert!(sparse_grid.sparse);
    let expected = brute_pairs(&positions, 4.0);
    assert_eq!(listed_pairs(&dense, positions.len(), 4.0), expected);
    assert_eq!(listed_pairs(&sparse, positions.len(), 4.0), expected);
    let probe = [30.0, 30.0, 30.0];
    let (mut a, mut b) = (Vec::new(), Vec::new());
    dense
        .for_each_neighbor(probe, 4.0, |t, _| a.push(t))
        .unwrap();
    sparse
        .for_each_neighbor(probe, 4.0, |t, _| b.push(t))
        .unwrap();
    a.sort_unstable();
    b.sort_unstable();
    assert_eq!(a, b);
}

#[test]
#[allow(clippy::cast_possible_truncation)]
fn a_huge_extent_index_fits_a_small_execution_budget() {
    use molframe_core::{ExecutionContext, MemoryBudget, ScratchPolicy};
    let mut positions = scatter(2000, 40.0, 3);
    positions.push([5.0e7, 5.0e7, 5.0e7]);
    let targets: Vec<u32> = (0..positions.len() as u32).collect();
    let context = ExecutionContext::builder()
        .memory_budget(MemoryBudget::new(1 << 20).unwrap())
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .unwrap();
    let list = CellList::build_in(
        &positions,
        &targets,
        3.0,
        None,
        CellGridOptions::default(),
        &context,
    )
    .expect("memory is proportional to atoms, not to extent");
    assert!(context.reserved_bytes() <= 1 << 20);
    assert_eq!(
        listed_pairs(&list, positions.len(), 3.0),
        brute_pairs(&positions, 3.0)
    );
}
