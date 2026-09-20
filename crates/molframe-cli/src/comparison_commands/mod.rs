//! Structure comparison command projections.

mod command;
mod render;

pub(super) use command::{ComparisonOptions, compare, map_chains, rmsd, superpose};

#[cfg(test)]
#[path = "comparison_commands_tests.rs"]
mod tests;
