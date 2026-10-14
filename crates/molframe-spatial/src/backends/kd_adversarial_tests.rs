//! Exactness of the bucketed k-d tree against brute force on hostile layouts.

use super::*;
use crate::brute::{distance_squared, pairs_same_selection};
use molframe_core::numeric::{f64_to_f32, u64_to_f64, usize_to_f64, usize_to_u32};
use molframe_core::structure::UnitCell;

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> f32 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        f64_to_f32(u64_to_f64(self.0 >> 40) / 16_777_216.0)
    }

    fn point(&mut self, extent: [f32; 3]) -> [f32; 3] {
        [
            self.next() * extent[0],
            self.next() * extent[1],
            self.next() * extent[2],
        ]
    }
}

fn layouts() -> Vec<(&'static str, Vec<[f32; 3]>, f32)> {
    let mut rng = Lcg(7);
    let uniform: Vec<_> = (0..3000).map(|_| rng.point([40.0; 3])).collect();
    let needle: Vec<_> = (0..2000).map(|_| rng.point([200.0, 0.02, 0.02])).collect();
    let slab: Vec<_> = (0..3000).map(|_| rng.point([60.0, 60.0, 0.4])).collect();
    let mut clustered = Vec::new();
    for centre in 0..6_i16 {
        let origin = [f32::from(centre) * 30.0, 0.0, 0.0];
        for _ in 0..500 {
            let local = rng.point([1.5; 3]);
            clustered.push([origin[0] + local[0], local[1], local[2]]);
        }
    }
    clustered.extend((0..400).map(|_| rng.point([300.0, 100.0, 100.0])));
    let offset: Vec<_> = uniform
        .iter()
        .map(|p| [p[0] + 1.0e5, p[1] + 1.0e5, p[2] + 1.0e5])
        .collect();
    let collinear: Vec<_> = (0..2500_u16)
        .map(|i| [f32::from(i) * 0.37, 1.0, 2.0])
        .collect();
    let mut outliers = uniform.clone();
    for (slot, sign) in [(0, 1.0), (1, -1.0), (2, 1.0)] {
        outliers[slot] = [sign * 1.0e6, 3.0 * sign, -1.0e6];
    }
    let duplicates: Vec<_> = (0..2400).map(|i| uniform[i % 600]).collect();
    vec![
        ("uniform", uniform, 3.0),
        ("needle", needle, 0.7),
        ("slab", slab, 2.0),
        ("clustered", clustered, 1.0),
        ("offset", offset, 3.0),
        ("collinear", collinear, 1.0),
        ("outliers", outliers, 3.0),
        ("duplicates", duplicates, 3.0),
    ]
}

fn triclinic() -> PeriodicBox {
    PeriodicBox::from_cell(UnitCell {
        lengths: [18.0, 21.0, 24.0],
        angles: [70.0, 80.0, 65.0],
    })
    .expect("cell")
}

fn all(count: usize) -> Vec<u32> {
    (0..usize_to_u32(count)).collect()
}

fn key(pair: &NeighborPair) -> (u32, u32, u32) {
    (pair.first, pair.second, pair.distance_squared.to_bits())
}

fn context(workers: usize) -> ExecutionContext {
    ExecutionContext::builder()
        .worker_budget(workers)
        .build()
        .expect("context")
}

fn check_pairs(name: &str, positions: &[[f32; 3]], cutoff: f32, periodic: Option<PeriodicBox>) {
    let atoms = all(positions.len());
    let expected = pairs_same_selection(positions, &atoms, cutoff * cutoff, periodic.as_ref());
    for workers in [1, 4] {
        let context = context(workers);
        let tree = KdTree::build_in(
            positions,
            &atoms,
            periodic,
            KdPeriodicOptions::default(),
            &context,
        )
        .expect("tree");
        let found = tree
            .pairs_in(&atoms, cutoff, true, true, &context)
            .expect("pairs");
        assert_eq!(
            found.iter().map(key).collect::<Vec<_>>(),
            expected.iter().map(key).collect::<Vec<_>>(),
            "{name} pairs with {workers} workers"
        );
        let direct = tree.pairs(&atoms, cutoff).expect("serial pairs");
        assert_eq!(direct.len(), expected.len(), "{name} serial pair count");
    }
}

#[test]
fn radius_pairs_equal_brute_force_on_adversarial_layouts() {
    for (name, positions, cutoff) in layouts() {
        check_pairs(name, &positions, cutoff, None);
    }
}

#[test]
fn periodic_triclinic_pairs_equal_brute_force() {
    let mut rng = Lcg(99);
    // Points deliberately outside the primary cell exercise wrapping.
    let positions: Vec<_> = (0..1500)
        .map(|_| {
            let p = rng.point([60.0, 60.0, 60.0]);
            [p[0] - 20.0, p[1] - 20.0, p[2] - 20.0]
        })
        .collect();
    check_pairs("triclinic", &positions, 4.0, Some(triclinic()));
}

fn check_knn(
    name: &str,
    positions: &[[f32; 3]],
    k: usize,
    max: Option<f32>,
    periodic: Option<PeriodicBox>,
) {
    let atoms = all(positions.len());
    let reference = |atom: u32| {
        let point = positions[atom as usize];
        let mut found: Vec<(f32, u32)> = atoms
            .iter()
            .filter(|&&t| t != atom)
            .map(|&t| {
                (
                    distance_squared(point, positions[t as usize], periodic.as_ref()),
                    t,
                )
            })
            .filter(|&(d, _)| max.is_none_or(|m| d <= m * m))
            .collect();
        found.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        found.truncate(k);
        found
    };
    let options = KnnOptions {
        max_distance: max,
        exclude_self: true,
    };
    let mut first: Option<KnnResult> = None;
    for workers in [1, 4] {
        let context = context(workers);
        let tree = KdTree::build_in(
            positions,
            &atoms,
            periodic,
            KdPeriodicOptions::default(),
            &context,
        )
        .expect("tree");
        let result = tree.knn(&atoms, k, options, &context).expect("knn");
        assert_eq!(result.k, k);
        for (row, &atom) in atoms.iter().enumerate().step_by(7) {
            let expected = reference(atom);
            let ids = &result.indices[row * k..(row + 1) * k];
            let d2 = &result.distances_squared[row * k..(row + 1) * k];
            for (slot, &(distance, id)) in expected.iter().enumerate() {
                assert_eq!(ids[slot], id, "{name} row {row} slot {slot}");
                assert_eq!(d2[slot].to_bits(), distance.to_bits(), "{name} distance");
            }
            for slot in expected.len()..k {
                assert_eq!(ids[slot], NO_NEIGHBOR, "{name} padding");
                assert!(d2[slot].is_infinite());
            }
        }
        if let Some(previous) = &first {
            assert_eq!(previous, &result, "{name} differs across worker counts");
        }
        first = Some(result);
    }
}

#[test]
fn knn_equals_brute_force_on_adversarial_layouts() {
    for (name, positions, cutoff) in layouts() {
        check_knn(name, &positions, 9, None, None);
        check_knn(name, &positions, 9, Some(cutoff), None);
    }
}

#[test]
fn periodic_knn_equals_brute_force_minimum_image() {
    let mut rng = Lcg(5);
    let positions: Vec<_> = (0..900)
        .map(|_| {
            let p = rng.point([40.0, 40.0, 40.0]);
            [p[0] - 10.0, p[1] - 10.0, p[2] - 10.0]
        })
        .collect();
    check_knn("triclinic", &positions, 8, None, Some(triclinic()));
    check_knn(
        "triclinic bounded",
        &positions,
        8,
        Some(5.0),
        Some(triclinic()),
    );
}

#[test]
fn knn_pads_short_rows_and_skips_non_finite_points() {
    let positions = [[0.0; 3], [1.0, 0.0, 0.0], [f32::NAN; 3]];
    let context = context(2);
    let atoms = all(3);
    let tree = KdTree::build_in(
        &positions,
        &atoms,
        None,
        KdPeriodicOptions::default(),
        &context,
    )
    .expect("tree");
    let options = KnnOptions {
        max_distance: None,
        exclude_self: true,
    };
    let result = tree.knn(&atoms, 3, options, &context).expect("knn");
    assert_eq!(result.indices[..3], [1, NO_NEIGHBOR, NO_NEIGHBOR]);
    assert_eq!(result.indices[3..6], [0, NO_NEIGHBOR, NO_NEIGHBOR]);
    assert_eq!(result.indices[6..], [NO_NEIGHBOR; 3]);
    assert!(result.distances_squared[0].total_cmp(&1.0).is_eq());
}

#[test]
fn tree_depth_stays_logarithmic_for_degenerate_inputs() {
    for (name, positions, _) in layouts() {
        let tree = KdTree::build(&positions, &all(positions.len()), None).expect("tree");
        let leaves = usize_to_f64(tree.split_value.len() + 1);
        assert!(
            leaves <= 2.0 * usize_to_f64(positions.len()) / 16.0 + 1.0,
            "{name}: {leaves} leaves for {}",
            positions.len()
        );
    }
}
