use super::{PoreProfileOptions, pore_profile};

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
    };
    assert!(pore_profile(&[[0.0; 3]], &[1.0], options).is_err());
}
