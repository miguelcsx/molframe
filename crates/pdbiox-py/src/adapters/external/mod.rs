//! External object-model exports backed by one shared Rust projection.

mod array_models;
mod common;
mod hierarchy_models;
mod molecular_dynamics;
mod openmm;
mod openstructure;
mod parmed;
mod rdkit;

#[cfg(test)]
#[path = "external_tests.rs"]
mod tests;
