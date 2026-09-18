use super::{AtomDepthOptions, atom_depths};

fn depths(atoms: &[[f32; 3]], surface: &[[f32; 3]]) -> Vec<f32> {
    match atom_depths(atoms, surface, AtomDepthOptions { cell_size: 3.0 }) {
        Ok(depths) => depths,
        Err(error) => panic!("explicit valid options failed: {error}"),
    }
}

#[test]
fn an_atom_takes_the_distance_to_its_nearest_surface_point() {
    let depths = depths(&[[0.0, 0.0, 0.0]], &[[0.0, 0.0, 5.0], [10.0, 0.0, 0.0]]);
    assert!((depths[0] - 5.0).abs() < 1e-4, "depth {}", depths[0]);
}

#[test]
fn an_atom_on_the_surface_has_zero_depth() {
    let depths = depths(&[[1.0, 2.0, 3.0]], &[[1.0, 2.0, 3.0]]);
    assert!(depths[0] < 1e-6);
}

#[test]
fn the_nearest_point_is_found_even_when_it_is_far() {
    // The nearest point is 20 Å away, several grid shells out; the search must
    // still reach it rather than report infinity.
    let depths = depths(&[[0.0, 0.0, 0.0]], &[[20.0, 0.0, 0.0]]);
    assert!((depths[0] - 20.0).abs() < 1e-3, "depth {}", depths[0]);
}

#[test]
fn without_a_surface_every_depth_is_infinite() {
    let depths = depths(&[[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]], &[]);
    assert!(depths.iter().all(|depth| depth.is_infinite()));
}
