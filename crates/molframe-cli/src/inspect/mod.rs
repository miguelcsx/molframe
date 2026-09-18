//! General, read-only inspection commands.

mod command;

pub(super) use command::{assemblies, categories, sequence, show};

#[cfg(test)]
#[path = "inspect_tests.rs"]
mod tests;
