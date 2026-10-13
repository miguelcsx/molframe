use super::{Grid, Point};

fn near(grid: &Grid, point: Point) -> Vec<u32> {
    let mut found = Vec::new();
    grid.for_each_near(point, |index, _| found.push(index));
    found.sort_unstable();
    found
}

fn brute(items: &[(u32, Point)], grid: &Grid, point: Point) -> Vec<u32> {
    let centre = grid.key(point);
    let mut found: Vec<u32> = items
        .iter()
        .filter(|(_, other)| {
            let key = grid.key(*other);
            (0..3).all(|axis| (key[axis] - centre[axis]).abs() <= 1)
        })
        .map(|(index, _)| *index)
        .collect();
    found.sort_unstable();
    found
}

fn scatter(count: u32, edge: f32) -> Vec<(u32, Point)> {
    let mut state = 7u64;
    let mut next = || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        f32::from(u16::try_from(state >> 48).unwrap_or(0)) / 65_536.0 * edge
    };
    (0..count).map(|i| (i, [next(), next(), next()])).collect()
}

#[test]
fn dense_neighbourhood_matches_brute_force() {
    let items = scatter(500, 40.0);
    let grid = Grid::build(4.0, items.iter().copied());
    for (_, point) in &items {
        assert_eq!(near(&grid, *point), brute(&items, &grid, *point));
    }
}

#[test]
fn sparse_far_flung_neighbourhood_matches_brute_force() {
    let mut items = scatter(200, 20.0);
    items.push((1000, [1.0e9, -1.0e9, 5.0e8]));
    items.push((1001, [1.0e9, -1.0e9, 5.0e8]));
    items.push((1002, [f32::MAX, f32::MIN, 0.0]));
    let grid = Grid::build(3.0, items.iter().copied());
    for (_, point) in &items {
        assert_eq!(near(&grid, *point), brute(&items, &grid, *point));
    }
    assert!(near(&grid, [-5.0e8, 0.0, 0.0]).is_empty());
}

#[test]
fn infinite_cell_collects_everything_and_empty_grid_finds_nothing() {
    let items = scatter(50, 100.0);
    let grid = Grid::build(f64::INFINITY, items.iter().copied());
    assert_eq!(near(&grid, [1.0, 1.0, 1.0]).len(), 50);
    let empty = Grid::build(2.0, std::iter::empty());
    assert!(near(&empty, [0.0; 3]).is_empty());
}
