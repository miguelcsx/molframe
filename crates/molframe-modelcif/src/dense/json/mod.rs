//! `AlphaFold` predicted-aligned-error JSON import into a dense matrix.

mod import;
mod scan;

pub use import::read_pae_json;
