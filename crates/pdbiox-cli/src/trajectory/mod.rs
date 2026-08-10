//! Declarative CLI projection of native trajectory I/O.

mod command;
mod summary;

pub(super) use command::{convert, extract, info, rmsd};

#[cfg(test)]
#[path = "trajectory_tests.rs"]
mod tests;
