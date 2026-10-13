use super::*;
use crate::fibonacci_sphere;
use molframe_core::{MemoryBudget, ScratchPolicy};

#[test]
fn streamed_areas_match_resident_bits_at_every_worker_count() {
    let positions: Vec<_> = (0_u16..150)
        .map(|index| {
            [
                f32::from(index % 11) * 1.7,
                f32::from(index / 11) * 1.3,
                f32::from(index % 3) * 0.4,
            ]
        })
        .collect();
    let radii: Vec<_> = (0_u16..150)
        .map(|index| 1.0 + f32::from(index % 5) * 0.1)
        .collect();
    let reference =
        crate::shrake_rupley(&positions, &radii, 1.4, 200, &ExecutionContext::default())
            .expect("reference");
    for workers in [1, 2, 4, 8] {
        let context = ExecutionContext::builder()
            .worker_budget(workers)
            .build()
            .expect("context");
        let mut result = Vec::new();
        visit_shrake_rupley(
            &positions,
            &radii,
            1.4,
            200,
            None,
            &context,
            |atom, area| {
                assert_eq!(atom, result.len());
                result.push(area);
                Ok(())
            },
        )
        .expect("stream");
        assert_eq!(
            result.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            reference.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
        );
        assert_eq!(context.reserved_bytes(), 0);
    }
}

#[test]
fn dense_neighbourhoods_do_not_retain_quadratic_adjacency() {
    let positions = vec![[0.0; 3]; 1024];
    let radii = vec![1.0; 1024];
    let context = ExecutionContext::builder()
        .worker_budget(8)
        .memory_budget(MemoryBudget::new(131_072).expect("budget"))
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("context");
    let mut count = 0;
    visit_shrake_rupley(&positions, &radii, 1.4, 32, None, &context, |_, _| {
        count += 1;
        Ok(())
    })
    .expect("stream");
    assert_eq!(count, 1024);
    assert!(context.peak_reserved_bytes() <= 131_072);
    assert_eq!(context.reserved_bytes(), 0);
}

#[test]
fn refusal_cancellation_and_sink_failure_release_the_complete_workspace() {
    let context = ExecutionContext::builder()
        .memory_budget(MemoryBudget::new(1).expect("budget"))
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("context");
    assert!(matches!(
        visit_shrake_rupley(&[[0.0; 3]], &[1.0], 1.4, 32, None, &context, |_, _| panic!(
            "no capacity"
        )),
        Err(SasaError::Spatial(SpatialError::Memory(_)))
    ));
    assert_eq!(context.reserved_bytes(), 0);
    let context = ExecutionContext::default();
    let mut count = 0;
    assert!(
        visit_shrake_rupley(
            &[[0.0; 3]; 256],
            &[1.0; 256],
            1.4,
            32,
            None,
            &context,
            |_, _| {
                count += 1;
                Err(SasaError::from(SpatialError::Cancelled))
            }
        )
        .is_err()
    );
    assert_eq!(count, 1);
    assert_eq!(context.reserved_bytes(), 0);
    context.cancellation().cancel();
    assert!(matches!(
        visit_shrake_rupley(&[[0.0; 3]], &[1.0], 1.4, 32, None, &context, |_, _| panic!(
            "cancelled"
        )),
        Err(SasaError::Spatial(SpatialError::Cancelled))
    ));
}

#[test]
fn tangent_and_zero_radius_atoms_preserve_the_reference() {
    let positions = [[0.0; 3], [2.0, 0.0, 0.0], [3.0, 5.0, 0.0], [4.0, 0.0, 0.0]];
    let radii = [1.0, 1.0, 1.0, 0.0];
    let context = ExecutionContext::default();
    let reference = crate::shrake_rupley(&positions, &radii, 0.0, 65, &context).expect("reference");
    let mut result = Vec::new();
    visit_shrake_rupley(&positions, &radii, 0.0, 65, None, &context, |_, area| {
        result.push(area);
        Ok(())
    })
    .expect("stream");
    assert_eq!(result, reference);
}

#[test]
fn nonfinite_positions_are_refused_rather_than_reported_exposed() {
    let positions = [[0.0; 3], [2.0, 0.0, 0.0], [f32::NAN; 3]];
    let radii = [1.0; 3];
    let context = ExecutionContext::default();
    let result = visit_shrake_rupley(&positions, &radii, 0.0, 65, None, &context, |_, _| {
        panic!("no area may be emitted for refused input")
    });
    assert!(matches!(
        result,
        Err(SasaError::NonFinitePosition { atom: 2 })
    ));
}

#[test]
fn periodic_sampling_matches_a_direct_minimum_image_reference() {
    use molframe_core::structure::UnitCell;
    for angles in [[90.0; 3], [90.0, 90.0, 60.0]] {
        let periodic = PeriodicBox::from_cell(UnitCell {
            lengths: [10.0; 3],
            angles,
        })
        .expect("box");
        let positions = [[0.2, 0.3, 0.4], [9.8, 0.3, 0.4], [4.0, 4.0, 4.0]];
        let radii = [1.0, 1.2, 1.3];
        let directions = fibonacci_sphere(200);
        let expected: Vec<_> = positions
            .iter()
            .enumerate()
            .map(|(atom, point)| {
                let radius = f64::from(radii[atom]) + 1.4;
                let mut hidden = 0;
                let mut overlap = false;
                for (other, centre) in positions.iter().enumerate() {
                    if atom == other {
                        continue;
                    }
                    let reach = radius + f64::from(radii[other]) + 1.4;
                    overlap |= squared_distance(
                        [0.0; 3],
                        periodic.displacement_f64(point.map(f64::from), centre.map(f64::from)),
                    ) < reach * reach;
                }
                for &direction in &directions {
                    let sample = point_on_sphere(point.map(f64::from), radius, direction);
                    let covered = positions.iter().enumerate().any(|(other, centre)| {
                        let other_radius = f64::from(radii[other]) + 1.4;
                        atom != other
                            && squared_distance(
                                [0.0; 3],
                                periodic.displacement_f64(sample, centre.map(f64::from)),
                            ) < other_radius * other_radius
                    });
                    hidden += u16::from(covered);
                }
                if overlap {
                    f64::from(200 - hidden)
                        * (4.0 * core::f64::consts::PI / 200.0)
                        * radius
                        * radius
                } else {
                    4.0 * core::f64::consts::PI * radius * radius
                }
            })
            .collect();
        let mut actual = Vec::new();
        visit_shrake_rupley(
            &positions,
            &radii,
            1.4,
            200,
            Some(&periodic),
            &ExecutionContext::default(),
            |_, area| {
                actual.push(area);
                Ok(())
            },
        )
        .expect("periodic stream");
        for (actual, expected) in actual.iter().zip(expected) {
            assert!((actual - expected).abs() < 1e-5);
        }
    }
}

#[test]
fn sampler_reuses_its_charge_and_collected_areas_keep_their_final_owner() {
    use std::sync::Arc;
    let context = ExecutionContext::default();
    let radii = [1.0; 3];
    let sampler = SasaSampler::new(&radii, 1.4, 64, &context).expect("sampler");
    let shared_bytes = 3 * 4 + 64 * 24;
    assert_eq!(context.reserved_bytes(), shared_bytes);
    for _ in 0..100 {
        sampler
            .visit(&[[0.0; 3]; 3], None, |_, _| Ok(()))
            .expect("visit");
        assert_eq!(context.reserved_bytes(), shared_bytes);
    }
    let areas = Arc::new(sampler.collect(&[[0.0; 3]; 3], None).expect("collect"));
    assert_eq!(context.reserved_bytes(), shared_bytes + 3 * 8);
    let retained = Arc::clone(&areas);
    drop(sampler);
    drop(areas);
    assert_eq!(context.reserved_bytes(), 3 * 8);
    assert_eq!(retained.len(), 3);
    drop(retained);
    assert_eq!(context.reserved_bytes(), 0);
}

#[allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss
)]
fn lattice(count: usize, spacing: f32) -> (Vec<[f32; 3]>, Vec<f32>) {
    let side = (count as f64).cbrt().ceil() as usize;
    let positions = (0..count)
        .map(|index| {
            [
                (index % side) as f32 * spacing,
                (index / side % side) as f32 * spacing,
                (index / (side * side)) as f32 * spacing,
            ]
        })
        .collect();
    (positions, vec![1.7; count])
}

fn explicit(bytes: usize) -> ExecutionContext {
    ExecutionContext::builder()
        .worker_budget(2)
        .memory_budget(MemoryBudget::new(bytes).expect("budget"))
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .expect("context")
}

fn bits(areas: &[f64]) -> Vec<u64> {
    areas.iter().map(|area| area.to_bits()).collect()
}

#[test]
fn a_budget_too_small_for_one_index_falls_back_to_slabs_with_identical_bits() {
    let (positions, radii) = lattice(6000, 3.0);
    let reference =
        crate::shrake_rupley(&positions, &radii, 1.4, 24, &explicit(1 << 30)).expect("reference");
    // Smallest budget at which the sweep still succeeds, by bisection.
    let (mut low, mut high) = (0_usize, 1_usize << 30);
    while high - low > 1024 {
        let middle = low + (high - low) / 2;
        if crate::shrake_rupley(&positions, &radii, 1.4, 24, &explicit(middle)).is_ok() {
            high = middle;
        } else {
            low = middle;
        }
    }
    // The whole-structure index needs far more than the slab sweep does.
    let whole = (positions.len() * 60).max(high + 1);
    assert!(
        high < whole / 2,
        "slab minimum {high} vs whole index {whole}"
    );
    let slabbed =
        crate::shrake_rupley(&positions, &radii, 1.4, 24, &explicit(high)).expect("slabbed");
    assert_eq!(bits(&slabbed), bits(&reference));
}

#[test]
fn an_explicit_tiny_budget_still_refuses_with_the_remedy() {
    let (positions, radii) = lattice(2000, 3.0);
    let error = crate::shrake_rupley(&positions, &radii, 1.4, 24, &explicit(2000))
        .expect_err("refused")
        .to_string();
    for needle in [
        "bytes",
        "memory_budget",
        "--memory-budget",
        "MOLFRAME_MEMORY_BUDGET",
    ] {
        assert!(error.contains(needle), "{needle}: {error}");
    }
}

#[test]
#[ignore = "3M atoms: run with `cargo test --release -p molframe-surface -- --ignored`"]
fn three_million_atoms_run_under_the_default_budget() {
    let (positions, radii) = lattice(3_000_000, 3.2);
    let context = ExecutionContext::default();
    let mut atoms = 0;
    let mut total = 0.0;
    visit_shrake_rupley(&positions, &radii, 1.4, 8, None, &context, |_, area| {
        atoms += 1;
        total += area;
        Ok(())
    })
    .expect("the default budget must not refuse");
    assert_eq!(atoms, 3_000_000);
    assert!(total > 0.0);
}
