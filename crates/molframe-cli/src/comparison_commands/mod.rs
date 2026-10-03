//! Structure comparison command projections.

mod ce;
mod command;
mod docking;
mod overlap;
mod regions;
mod render;

pub(super) use command::{ComparisonOptions, MappedScoring, compare, map_chains, rmsd, superpose};

#[cfg(test)]
#[path = "comparison_commands_tests.rs"]
mod tests;
