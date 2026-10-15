//! Commands for complete topology-and-coordinate system containers.

mod command;

pub(crate) use command::{copy, info};

#[cfg(test)]
#[path = "system_tests.rs"]
mod tests;
