use super::{SasaError, shrake_rupley, surface_points};
use core::f64::consts::PI;
use molframe_core::ExecutionContext;

fn context() -> ExecutionContext {
    ExecutionContext::default()
}

fn sphere_area(radius: f64) -> f64 {
    4.0 * PI * radius * radius
}

// The kernel widens each radius with the probe in f64 after an f32 input, so a
// reference area must expand the same way rather than assume exact decimals.
fn expanded(radius: f32, probe: f32) -> f64 {
    f64::from(radius) + f64::from(probe)
}

#[test]
fn a_lone_atom_exposes_its_whole_expanded_sphere() {
    let Ok(areas) = shrake_rupley(&[[0.0, 0.0, 0.0]], &[1.5], 1.4, 400, &context()) else {
        panic!("a single atom is valid");
    };
    assert!((areas[0] - sphere_area(expanded(1.5, 1.4))).abs() < 1e-6);
}

#[test]
fn a_small_atom_buried_inside_a_larger_one_has_no_accessible_area() {
    // Concentric: the small sphere lies wholly within the large one.
    let positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
    let Ok(areas) = shrake_rupley(&positions, &[1.0, 3.0], 0.0, 400, &context()) else {
        panic!("valid");
    };
    assert!(areas[0].abs() < 1e-9, "buried atom kept area {}", areas[0]);
    assert!((areas[1] - sphere_area(3.0)).abs() < 1e-6);
}

#[test]
fn atoms_far_apart_do_not_shade_each_other() {
    let positions = [[0.0, 0.0, 0.0], [100.0, 0.0, 0.0]];
    let Ok(areas) = shrake_rupley(&positions, &[1.5, 1.5], 1.4, 400, &context()) else {
        panic!("valid");
    };
    let full = sphere_area(expanded(1.5, 1.4));
    assert!((areas[0] - full).abs() < 1e-6);
    assert!((areas[1] - full).abs() < 1e-6);
}

#[test]
fn touching_atoms_lose_area_and_share_the_loss() {
    // Two equal spheres whose expanded radii overlap occlude each other equally.
    // The physical areas are identical; a finite, non-mirror-symmetric point
    // lattice only reproduces that to within its sampling, so the two match
    // closely rather than exactly, and both sit below a full sphere.
    let positions = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let Ok(areas) = shrake_rupley(&positions, &[1.5, 1.5], 0.0, 2000, &context()) else {
        panic!("valid");
    };
    let full = sphere_area(1.5);
    assert!(
        (areas[0] - areas[1]).abs() < 0.01 * full,
        "occlusion should be near-symmetric: {} vs {}",
        areas[0],
        areas[1]
    );
    assert!(areas[0] < full && areas[0] > 0.0, "expected partial area");
    assert!(areas[1] < full && areas[1] > 0.0, "expected partial area");
}

#[test]
fn a_mismatched_radius_count_is_rejected() {
    let Err(error) = shrake_rupley(&[[0.0, 0.0, 0.0]], &[1.0, 2.0], 1.4, 100, &context()) else {
        panic!("mismatch should be rejected");
    };
    assert!(matches!(error, SasaError::LengthMismatch { .. }));
}

#[test]
fn a_negative_probe_is_rejected() {
    let Err(error) = shrake_rupley(&[[0.0, 0.0, 0.0]], &[1.0], -1.0, 100, &context()) else {
        panic!("negative probe should be rejected");
    };
    assert!(matches!(error, SasaError::InvalidProbe));
}

#[test]
fn zero_test_points_is_rejected() {
    let Err(error) = shrake_rupley(&[[0.0, 0.0, 0.0]], &[1.0], 1.4, 0, &context()) else {
        panic!("zero points should be rejected");
    };
    assert!(matches!(error, SasaError::NoPoints));
}

#[test]
fn an_empty_set_has_no_areas() {
    let Ok(areas) = shrake_rupley(&[], &[], 1.4, 100, &context()) else {
        panic!("an empty set is valid");
    };
    assert!(areas.is_empty());
}

#[test]
fn a_lone_atom_yields_one_surface_point_per_sample() {
    let Ok(points) = surface_points(&[[0.0, 0.0, 0.0]], &[1.5], 1.4, 200, &context()) else {
        panic!("valid");
    };
    assert_eq!(points.len(), 200);
    for point in &points {
        assert_eq!(point.atom, 0);
        // Every point sits at the expanded radius from the centre.
        let radius = expanded(1.5, 1.4);
        let distance = (f64::from(point.position[0]).powi(2)
            + f64::from(point.position[1]).powi(2)
            + f64::from(point.position[2]).powi(2))
        .sqrt();
        assert!((distance - radius).abs() < 1e-4, "distance {distance}");
    }
}

#[test]
fn a_buried_atom_contributes_no_surface_points() {
    // A small atom concentric inside a large one is fully buried.
    let positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
    let Ok(points) = surface_points(&positions, &[1.0, 3.0], 0.0, 400, &context()) else {
        panic!("valid");
    };
    assert!(
        points.iter().all(|point| point.atom != 0),
        "atom 0 is buried"
    );
}

// A packed grid so most atoms have neighbours and the parallel partition
// crosses buried, exposed and boundary atoms.
fn packed_grid() -> (Vec<[f32; 3]>, Vec<f32>) {
    let mut positions = Vec::new();
    for x in 0_u8..7 {
        for y in 0_u8..7 {
            for z in 0_u8..7 {
                positions.push([f32::from(x) * 1.6, f32::from(y) * 1.6, f32::from(z) * 1.6]);
            }
        }
    }
    let radii = vec![1.7; positions.len()];
    (positions, radii)
}

#[test]
fn shrake_rupley_is_bit_identical_across_worker_counts() {
    let (positions, radii) = packed_grid();
    let Ok(serial) = shrake_rupley(&positions, &radii, 1.4, 200, &context()) else {
        panic!("serial run is valid");
    };
    for workers in [2, 3, 4, 8, 16] {
        let worker_context = ExecutionContext::builder()
            .worker_budget(workers)
            .build()
            .expect("worker context is valid");
        let Ok(parallel) = shrake_rupley(&positions, &radii, 1.4, 200, &worker_context) else {
            panic!("parallel run with {workers} workers is valid");
        };
        assert_eq!(
            serial.to_bits_vec(),
            parallel.to_bits_vec(),
            "worker count {workers} changed the result"
        );
    }
}

#[test]
fn surface_points_are_identical_across_worker_counts() {
    let (positions, radii) = packed_grid();
    let Ok(serial) = surface_points(&positions, &radii, 1.4, 96, &context()) else {
        panic!("serial run is valid");
    };
    for workers in [2, 4, 8, 16] {
        let worker_context = ExecutionContext::builder()
            .worker_budget(workers)
            .build()
            .expect("worker context is valid");
        let Ok(parallel) = surface_points(&positions, &radii, 1.4, 96, &worker_context) else {
            panic!("parallel run with {workers} workers is valid");
        };
        assert_eq!(
            serial, parallel,
            "worker count {workers} changed the points"
        );
    }
}

// Exact bit comparison so a reordered floating-point sum would fail, not round away.
trait ToBitsVec {
    fn to_bits_vec(&self) -> Vec<u64>;
}

impl ToBitsVec for Vec<f64> {
    fn to_bits_vec(&self) -> Vec<u64> {
        self.iter().map(|value| value.to_bits()).collect()
    }
}

fn budgeted(workers: usize, bytes: usize) -> ExecutionContext {
    ExecutionContext::builder()
        .worker_budget(workers)
        .memory_budget(molframe_core::MemoryBudget::new(bytes).expect("a positive budget"))
        .scratch_policy(molframe_core::ScratchPolicy::new(0))
        .build()
        .expect("a valid context")
}

#[test]
fn surface_points_stream_through_a_bounded_window_and_match_at_every_worker_count() {
    let (positions, radii) = packed_grid();
    let Ok(serial) = surface_points(&positions, &radii, 1.4, 64, &context()) else {
        panic!("serial run is valid");
    };
    assert!(!serial.is_empty());
    for workers in [1, 2, 4, 8] {
        // Room for the result and a few blocks in flight, not for every block at once.
        let roomy = budgeted(workers, 64 * 1024 * 1024);
        let Ok(points) = surface_points(&positions, &radii, 1.4, 64, &roomy) else {
            panic!("a roomy budget suffices at {workers} workers");
        };
        assert_eq!(points, serial, "{workers} workers");
        assert_eq!(
            roomy.reserved_bytes(),
            0,
            "{workers} workers release their reservation"
        );
    }
}

#[test]
fn a_budget_below_one_block_is_refused_instead_of_overcommitted() {
    let (positions, radii) = packed_grid();
    let tiny = budgeted(4, 1024);
    assert!(matches!(
        surface_points(&positions, &radii, 1.4, 64, &tiny),
        Err(SasaError::Spatial(molframe_spatial::SpatialError::Memory(
            _
        )))
    ));
}

mod non_finite_and_parallel {
    use super::{SasaError, context, shrake_rupley, surface_points};
    use crate::lee_richards::lee_richards;
    use crate::surface_points_at_density;
    use molframe_core::ExecutionContext;

    const BAD: [f32; 3] = [f32::NAN, 0.0, 0.0];

    fn cloud() -> (Vec<[f32; 3]>, Vec<f32>) {
        let positions: Vec<[f32; 3]> = (0..300_u16)
            .map(|index| {
                let step = f32::from(index);
                [
                    (step * 0.37).sin() * 9.0,
                    (step * 0.91).cos() * 9.0,
                    (step * 0.13).sin() * 9.0,
                ]
            })
            .collect();
        let radii = (0..300_u16).map(|i| 1.2 + f32::from(i % 5) * 0.2).collect();
        (positions, radii)
    }

    fn rejected(error: &SasaError, atom: usize) {
        assert!(
            matches!(error, SasaError::NonFinitePosition { atom: found } if *found == atom),
            "unexpected error {error:?}"
        );
    }

    #[test]
    fn every_kernel_refuses_non_finite_coordinates() {
        let ctx = context();
        for bad in [
            BAD,
            [0.0, f32::INFINITY, 0.0],
            [0.0, 0.0, f32::NEG_INFINITY],
        ] {
            let positions = [[0.0, 0.0, 0.0], bad];
            let radii = [1.7, 1.7];
            rejected(
                &shrake_rupley(&positions, &radii, 1.4, 100, &ctx).unwrap_err(),
                1,
            );
            rejected(
                &lee_richards(&positions, &radii, 1.4, 20, &ctx).unwrap_err(),
                1,
            );
            rejected(
                &surface_points(&positions, &radii, 1.4, 50, &ctx).unwrap_err(),
                1,
            );
            rejected(
                &surface_points_at_density(&positions, &radii, 1.4, 1.0, &ctx).unwrap_err(),
                1,
            );
        }
    }

    #[test]
    fn a_lone_non_finite_atom_is_not_reported_exposed() {
        let result = shrake_rupley(&[BAD], &[1.7], 1.4, 100, &context());
        assert!(matches!(
            result,
            Err(SasaError::NonFinitePosition { atom: 0 })
        ));
    }

    #[test]
    fn density_sampling_is_identical_at_any_worker_count() {
        let (positions, radii) = cloud();
        let run = |workers| {
            let ctx = ExecutionContext::builder()
                .worker_budget(workers)
                .build()
                .unwrap();
            surface_points_at_density(&positions, &radii, 1.4, 1.5, &ctx).unwrap()
        };
        let serial = run(1);
        assert!(!serial.is_empty());
        assert_eq!(serial, run(8));
        assert!(serial.windows(2).all(|pair| pair[0].atom <= pair[1].atom));
    }

    #[test]
    fn lee_richards_is_identical_at_any_worker_count_and_tracks_shrake_rupley() {
        let (positions, radii) = cloud();
        let run = |workers| {
            let ctx = ExecutionContext::builder()
                .worker_budget(workers)
                .build()
                .unwrap();
            lee_richards(&positions, &radii, 1.4, 60, &ctx).unwrap()
        };
        let areas = run(1);
        assert_eq!(areas, run(8));
        let reference = shrake_rupley(&positions, &radii, 1.4, 4000, &context()).unwrap();
        let (lee, sphere): (f64, f64) = (areas.iter().sum(), reference.iter().sum());
        assert!((lee - sphere).abs() / sphere < 0.02, "{lee} vs {sphere}");
    }
}
