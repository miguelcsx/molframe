//! Deterministic coordinate and radius fixture transformations.

use molframe_core::index::ModelIndex;
use molframe_core::structure::Structure;

/// The first model's atom positions as an owned slice.
#[must_use]
pub fn coordinates(structure: &Structure) -> Vec<[f32; 3]> {
    structure.positions().to_vec()
}

/// One model's atom positions as an owned slice, if the model exists.
#[must_use]
pub fn model_coordinates(structure: &Structure, model: usize) -> Option<Vec<[f32; 3]>> {
    let index = match u32::try_from(model) {
        Ok(index) => ModelIndex::new(index),
        Err(_) => return None,
    };
    structure.model_positions(index).map(<[[f32; 3]]>::to_vec)
}

/// A deterministic small perturbation of a coordinate set.
#[must_use]
pub fn perturbed(coordinates: &[[f32; 3]], jitter: f32) -> Vec<[f32; 3]> {
    coordinates
        .iter()
        .map(|&[x, y, z]| {
            [
                x + jitter * (x - x.floor() - 0.5),
                y + jitter * (y - y.floor() - 0.5),
                z + jitter * (z - z.floor() - 0.5),
            ]
        })
        .collect()
}

/// A uniform per-atom radius vector of the given length.
#[must_use]
pub fn uniform_radii(len: usize, radius: f64) -> Vec<f64> {
    vec![radius; len]
}
