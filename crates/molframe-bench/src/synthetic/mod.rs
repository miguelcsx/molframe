//! Generated fixtures for measurements larger than any committed structure.
//!
//! Scale is reached by tiling a real structure from a seed, so nothing large is
//! committed and nothing is downloaded. The byte sources emit format text as it
//! is consumed, which is what makes a measurement larger than memory possible.

mod emit;
mod format;
mod seed;
mod tile;

pub use emit::SyntheticCifSource;
pub use seed::Seed;
pub use tile::{Placement, Tile, TileAtom};
