//! Declarative CLI projection of native trajectory I/O.

mod command;
mod contacts;
mod rmsf;
mod streaming;
mod summary;
mod surface;

pub(crate) use command::{convert, extract, info, rmsd};
pub(crate) use contacts::contacts;
pub(crate) use rmsf::rmsf;
pub(crate) use surface::sasa;

#[cfg(test)]
#[path = "trajectory_tests.rs"]
mod tests;
