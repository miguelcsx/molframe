//! Commands for complete topology-and-coordinate system containers.

mod command;

pub(super) use command::{copy, info};

#[cfg(test)]
#[path = "system_tests.rs"]
mod tests;
