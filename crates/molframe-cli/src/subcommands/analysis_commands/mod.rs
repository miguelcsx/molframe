//! General analysis command projections.

mod command;

pub(crate) use command::{
    SasaWords, SseOptions, contacts, interfaces, neighbors, sasa, sse, sse_standard,
};

#[cfg(test)]
#[path = "analysis_commands_tests.rs"]
mod tests;
