//! General analysis command projections.

mod command;

pub(super) use command::{SseOptions, contacts, interfaces, neighbors, sasa, sse};

#[cfg(test)]
#[path = "analysis_commands_tests.rs"]
mod tests;
