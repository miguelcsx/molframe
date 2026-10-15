//! General, read-only inspection commands.

mod command;

pub(crate) use command::{assemblies, categories, sequence, show};

#[cfg(test)]
#[path = "inspect_tests.rs"]
mod tests;
