use super::{
    PoreError, PoreProfileOptions, PoreSample, add, cross, minimum_clearance, normalize,
    perpendicular, pore_profile, scale, transverse_coordinate,
};
use crate::numeric::{f32_to_usize, usize_to_f32};
use num_traits::ToPrimitive;

const MEMORY_LIMIT: usize = 100_000_000;

#[test]
fn the_memory_limit_defaults_to_one_hundred_megabytes_and_has_no_maximum() {
    let options = PoreProfileOptions::new(([0.0; 3], [0.0, 0.0, 1.0]), 0.0, 0.0, 1, 1.0, 0.25, 0.0);
    assert_eq!(
        options.memory_limit_bytes,
        PoreProfileOptions::DEFAULT_MEMORY_LIMIT_BYTES
    );
    assert_eq!(PoreProfileOptions::DEFAULT_MEMORY_LIMIT_BYTES, MEMORY_LIMIT);
    assert!(
        !matches!(
            pore_profile(
                &[[0.0; 3]],
                &[1.0],
                options.with_memory_limit(MEMORY_LIMIT * 100),
            ),
            Err(PoreError::InvalidOptions)
        ),
        "a caller who has provisioned the machine sets the ceiling, not the library"
    );
}

#[test]
fn symmetric_atoms_put_the_pore_centre_on_the_axis() {
    let positions = [
        [-3.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [0.0, -3.0, 0.0],
        [0.0, 3.0, 0.0],
    ];
    let options = PoreProfileOptions {
        axis_origin: [0.0; 3],
        axis_direction: [0.0, 0.0, 1.0],
        start: 0.0,
        end: 0.0,
        samples: 1,
        search_radius: 1.0,
        grid_spacing: 0.25,
        probe_radius: 0.0,
        memory_limit_bytes: MEMORY_LIMIT,
    };
    let Ok(profile) = pore_profile(&positions, &[1.0; 4], options) else {
        panic!("valid pore geometry");
    };
    assert!(
        profile[0]
            .centre
            .into_iter()
            .all(|value| value.abs() < f32::EPSILON)
    );
    assert!((profile[0].radius - 2.0).abs() < 1e-6);
}

#[test]
fn probe_radius_reduces_clearance_explicitly() {
    let options = PoreProfileOptions {
        axis_origin: [0.0; 3],
        axis_direction: [0.0, 0.0, 1.0],
        start: 0.0,
        end: 0.0,
        samples: 1,
        search_radius: 0.5,
        grid_spacing: 0.5,
        probe_radius: 0.4,
        memory_limit_bytes: MEMORY_LIMIT,
    };
    let positions = [
        [-2.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, -2.0, 0.0],
        [0.0, 2.0, 0.0],
    ];
    let Ok(profile) = pore_profile(&positions, &[1.0; 4], options) else {
        panic!("valid probe profile");
    };
    assert!((profile[0].radius - 0.6).abs() < 1e-6);
}

#[test]
fn no_axis_or_implicit_sampling_default_is_accepted() {
    let options = PoreProfileOptions {
        axis_origin: [0.0; 3],
        axis_direction: [0.0; 3],
        start: 0.0,
        end: 1.0,
        samples: 0,
        search_radius: 1.0,
        grid_spacing: 0.1,
        probe_radius: 0.0,
        memory_limit_bytes: MEMORY_LIMIT,
    };
    assert!(pore_profile(&[[0.0; 3]], &[1.0], options).is_err());
}

#[test]
fn spatial_pruning_is_exactly_equivalent_to_the_exhaustive_grid() {
    let mut positions = Vec::new();
    let mut radii = Vec::new();
    let mut state = 0x9e37_79b9_u32;
    for atom in 0..257 {
        state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        let x = unit(state) * 18.0 - 9.0;
        state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        let y = unit(state) * 18.0 - 9.0;
        state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        let z = unit(state) * 18.0 - 9.0;
        positions.push([x, y, z]);
        radii.push(0.8 + unit(state.rotate_left(atom % 31)));
    }
    let options = PoreProfileOptions {
        axis_origin: [0.25, -0.5, 0.75],
        axis_direction: [1.0, 2.0, -3.0],
        start: -3.0,
        end: 4.0,
        samples: 7,
        search_radius: 2.5,
        grid_spacing: 0.37,
        probe_radius: 1.1,
        memory_limit_bytes: MEMORY_LIMIT,
    };
    let Ok(actual) = pore_profile(&positions, &radii, options) else {
        panic!("indexed pore profile should fit its declared memory ceiling");
    };
    let expected = exhaustive_profile(&positions, &radii, options);
    assert_eq!(actual, expected);
}

#[test]
fn indexed_search_rejects_a_workspace_over_the_explicit_ceiling() {
    let positions = vec![[0.0, 0.0, 4.0]; 1_024];
    let options = PoreProfileOptions {
        axis_origin: [0.0; 3],
        axis_direction: [0.0, 0.0, 1.0],
        start: -2.0,
        end: 2.0,
        samples: 8,
        search_radius: 2.0,
        grid_spacing: 0.25,
        probe_radius: 0.0,
        memory_limit_bytes: 1_024,
    };
    assert!(matches!(
        pore_profile(&positions, &vec![1.0; positions.len()], options),
        Err(PoreError::MemoryLimit { required, limit: 1_024 }) if required > 1_024
    ));
}

#[test]
fn rounded_boundaries_and_first_maximum_ties_keep_exhaustive_order() {
    let positions = vec![[0.0; 3]; 64];
    let radii = vec![1.0; positions.len()];
    let options = PoreProfileOptions::new(([0.0; 3], [0.0, 0.0, 1.0]), 0.0, 0.0, 3, 1.0, 0.05, 0.0);
    let Ok(actual) = pore_profile(&positions, &radii, options) else {
        panic!("tie fixture should build its spatial index");
    };
    assert_eq!(actual, exhaustive_profile(&positions, &radii, options));
}

#[test]
fn extreme_finite_radii_do_not_invalidate_spatial_lower_bounds() {
    let mut positions = Vec::new();
    let mut radii = Vec::new();
    let radius_cycle = [0.0, f32::MIN_POSITIVE, f32::EPSILON, 1.0, f32::MAX / 8.0];
    for atom in 0_u16..33 {
        let offset = f32::from(atom % 11) - 5.0;
        positions.push([offset, f32::from(atom % 3) - 1.0, f32::from(atom % 5)]);
        radii.push(radius_cycle[usize::from(atom) % radius_cycle.len()]);
    }
    let options = PoreProfileOptions::new(([0.0; 3], [0.0, 0.0, 1.0]), 0.0, 1.0, 2, 1.0, 0.05, 0.0);
    let Ok(actual) = pore_profile(&positions, &radii, options) else {
        panic!("extreme finite radii should remain valid");
    };
    assert_eq!(actual, exhaustive_profile(&positions, &radii, options));
}

#[test]
fn derived_coordinate_overflow_matches_exhaustive_nan_handling() {
    let positions = vec![[0.0; 3]; 64];
    let radii = vec![1.0; positions.len()];
    let options = PoreProfileOptions::new(
        ([0.0; 3], [0.0, 0.0, 1.0]),
        -f32::MAX,
        f32::MAX,
        2,
        1.0,
        0.05,
        0.0,
    );
    let Ok(actual) = pore_profile(&positions, &radii, options) else {
        panic!("derived overflow remains representable by the exhaustive contract");
    };
    let expected = exhaustive_profile(&positions, &radii, options);
    assert_profile_bits_equal(&actual, &expected);
}

fn exhaustive_profile(
    positions: &[[f32; 3]],
    radii: &[f32],
    options: PoreProfileOptions,
) -> Vec<PoreSample> {
    let axis = normalize(options.axis_direction);
    let first_basis = perpendicular(axis);
    let second_basis = cross(axis, first_basis);
    let Some(steps) = f32_to_usize((2.0 * options.search_radius / options.grid_spacing).ceil())
    else {
        panic!("finite test grid");
    };
    let mut profile = Vec::with_capacity(options.samples);
    for sample in 0..options.samples {
        let fraction = if options.samples == 1 {
            0.0
        } else {
            usize_to_f32(sample) / usize_to_f32(options.samples - 1)
        };
        let axial_coordinate = options.start + fraction * (options.end - options.start);
        let axis_centre = add(options.axis_origin, scale(axis, axial_coordinate));
        let mut best_centre = axis_centre;
        let mut best_clearance = f32::NEG_INFINITY;
        for first_index in 0..=steps {
            let first = transverse_coordinate(first_index, steps, options.search_radius);
            for second_index in 0..=steps {
                let second = transverse_coordinate(second_index, steps, options.search_radius);
                if first.mul_add(first, second * second) > options.search_radius.powi(2) {
                    continue;
                }
                let candidate = add(
                    add(axis_centre, scale(first_basis, first)),
                    scale(second_basis, second),
                );
                let clearance =
                    minimum_clearance(candidate, positions, radii) - options.probe_radius;
                if clearance > best_clearance {
                    best_clearance = clearance;
                    best_centre = candidate;
                }
            }
        }
        profile.push(PoreSample {
            axial_coordinate,
            centre: best_centre,
            radius: best_clearance.max(0.0),
        });
    }
    profile
}

fn unit(value: u32) -> f32 {
    match (value >> 8).to_f32() {
        Some(value) => value / 16_777_215.0,
        None => 0.0,
    }
}

fn assert_profile_bits_equal(actual: &[PoreSample], expected: &[PoreSample]) {
    assert_eq!(actual.len(), expected.len());
    for (actual, expected) in actual.iter().zip(expected) {
        assert_eq!(
            actual.axial_coordinate.to_bits(),
            expected.axial_coordinate.to_bits()
        );
        assert_eq!(actual.radius.to_bits(), expected.radius.to_bits());
        for axis in 0..3 {
            assert_eq!(
                actual.centre[axis].to_bits(),
                expected.centre[axis].to_bits()
            );
        }
    }
}
