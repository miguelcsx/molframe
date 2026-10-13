//! Declarative scale-sweep infrastructure: size tiers, deterministic inputs
//! and the statistics that turn timings into a scaling exponent.

mod cloud;
mod rng;
mod stats;
mod tier;

pub use cloud::{Cloud, Ensemble, chain, ensemble, rotated, uniform_cloud};
pub use rng::SplitMix64;
pub use stats::{Checksum, Fit, fit_exponent};
pub use tier::{Tier, TierParseError};

#[cfg(test)]
#[path = "scale_tests.rs"]
mod tests;
