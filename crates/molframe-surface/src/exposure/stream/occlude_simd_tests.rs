use super::*;
use crate::neighbourhood::point_inside_sphere;

/// The scalar reference: a point is hidden by the first sphere that holds it.
fn scalar(points: &[[f64; 3]], neighbours: &[([f64; 3], f64)]) -> (Vec<u64>, u16) {
    let mut covered = vec![0_u64; points.len().div_ceil(64)];
    let mut hidden = 0;
    for (sample, point) in points.iter().enumerate() {
        if neighbours
            .iter()
            .any(|(centre, radius_squared)| point_inside_sphere(*point, *centre, *radius_squared))
        {
            covered[sample / 64] |= 1 << (sample % 64);
            hidden += 1;
        }
    }
    (covered, hidden)
}

fn lattice(count: usize, scale: f64) -> Vec<[f64; 3]> {
    (0..count)
        .map(|index| {
            let step = f64::from(u32::try_from(index).expect("small"));
            [
                (step * 0.37).sin() * scale,
                (step * 0.91).cos() * scale,
                (step * 0.13).sin() * scale,
            ]
        })
        .collect()
}

#[test]
fn the_tile_hides_exactly_the_points_the_scalar_test_hides() {
    let mut points = lattice(150, 3.0);
    points.push([f64::NAN, 0.0, 0.0]);
    points.push([0.0, f64::INFINITY, 0.0]);
    for length in [1, 3, 4, 63, 64] {
        let neighbours: Vec<([f64; 3], f64)> = lattice(length, 2.5)
            .into_iter()
            .enumerate()
            .map(|(index, centre)| {
                let radius = 0.5 + f64::from(u32::try_from(index % 7).expect("small")) * 0.3;
                (centre, radius * radius)
            })
            .collect();
        let mut tile = NeighbourTile::<64>::new();
        for (centre, radius_squared) in &neighbours {
            tile.push(*centre, *radius_squared);
        }
        let mut covered = vec![0_u64; points.len().div_ceil(64)];
        let hidden = tile.occlude(&points, &mut covered);
        assert_eq!(
            (covered, hidden),
            scalar(&points, &neighbours),
            "{length} spheres"
        );
    }
}

#[test]
fn a_cleared_tile_hides_nothing_and_previously_hidden_points_stay_counted_once() {
    let points = lattice(80, 1.0);
    let mut tile = NeighbourTile::<64>::new();
    tile.push([0.0; 3], 100.0);
    let mut covered = vec![0_u64; 2];
    assert_eq!(tile.occlude(&points, &mut covered), 80);
    assert_eq!(tile.occlude(&points, &mut covered), 0);
    tile.clear();
    tile.push([50.0; 3], 1.0);
    let mut fresh = vec![0_u64; 2];
    assert_eq!(tile.occlude(&points, &mut fresh), 0);
}
