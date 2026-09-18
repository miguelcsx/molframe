//! Declarative CLI projection of native trajectory I/O.

mod command;
mod contacts;
mod streaming;
mod summary;
mod surface;

pub(super) use command::{convert, extract, info, rmsd};
pub(super) use contacts::contacts;
pub(super) use surface::sasa;

#[cfg(test)]
#[path = "trajectory_tests.rs"]
mod tests;
