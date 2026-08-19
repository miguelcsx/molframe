//! Shared real-structure fixtures for the pdbiox benchmark suite.

mod coordinates;
mod fixtures;
mod samples;

pub use coordinates::{coordinates, model_coordinates, perturbed, uniform_radii};
pub use fixtures::{input, input_gzip, structure, structure_from_cif, structure_from_pdb};
pub use samples::{Sample, ccd_atp, ccd_hem, large_cif_gz};

#[cfg(test)]
#[path = "lib_tests.rs"]
mod tests;
