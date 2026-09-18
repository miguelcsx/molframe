//! Shared real-structure fixtures for the molframe benchmark suite.

mod coordinates;
mod fixtures;
mod numeric;
mod samples;
mod synthetic;

pub use coordinates::{coordinates, model_coordinates, perturbed, uniform_radii};
pub use fixtures::{input, input_gzip, structure, structure_from_cif, structure_from_pdb};
pub use samples::{Sample, ccd_atp, ccd_hem, large_cif_gz};
pub use synthetic::{Placement, Seed, SyntheticCifSource, Tile, TileAtom};

#[cfg(test)]
#[path = "lib_tests.rs"]
mod tests;
