use super::*;

/// Every cell's neighbourhood by exhaustive key lookup, the way the
/// cursor walk must reproduce it.
fn reference(grid: &CellGrid<u32>) -> Vec<(Range<usize>, Vec<Range<usize>>)> {
    let Layout::Sparse { cells } = &grid.layout else {
        panic!("the test builds sparse grids");
    };
    cells
        .iter()
        .map(|(base, own)| {
            let mut around = Vec::new();
            for dx in -1..=1_i32 {
                for dy in -1..=1_i32 {
                    for dz in -1..=1_i32 {
                        let key = CellKey([
                            base.0[0].saturating_add(dx),
                            base.0[1].saturating_add(dy),
                            base.0[2].saturating_add(dz),
                        ]);
                        if let Some((_, range)) = cells.iter().find(|(found, _)| *found == key) {
                            around.push(range.clone());
                        }
                    }
                }
            }
            (own.clone(), around)
        })
        .collect()
}

#[test]
fn sparse_neighbour_cursors_match_exhaustive_lookup() {
    let mut state = 7_u32;
    let mut next = move |span: i32| {
        state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        i32::try_from(state >> 20).expect("twelve bits") % span - span / 2
    };
    let mut entries = Vec::new();
    for item in 0..600_u32 {
        entries.push((CellKey([next(9), next(9), next(9)]), item));
    }
    // Far-apart cells force the sparse layout; extremes hit the fallback.
    for (item, key) in [
        [i32::MIN, 0, 0],
        [i32::MIN + 1, 0, 0],
        [i32::MAX, 5, -5],
        [i32::MAX - 1, 5, -5],
        [1_000_000_000, 1_000_000_000, 1_000_000_000],
        [-1_000_000_000, 3, 3],
    ]
    .into_iter()
    .enumerate()
    {
        entries.push((CellKey(key), 1_000 + u32::try_from(item).expect("small")));
    }
    let grid = CellGrid::sparse(entries);
    let expected = reference(&grid);
    let mut seen = Vec::new();
    grid.for_each_cell(|own, around| seen.push((own, around.to_vec())));
    assert_eq!(seen, expected);
    // Block-wise visits reproduce the whole sequence.
    let mut pieces = Vec::new();
    let count = grid.cell_count();
    for start in (0..count).step_by(37) {
        grid.visit_cells(start..start + 37, |own, around| {
            pieces.push((own, around.to_vec()));
        });
    }
    assert_eq!(pieces, expected);
}
