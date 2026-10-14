//! End-to-end CLI tests, built as one binary so the shared helpers in
//! `common` are compiled once and every helper is used.

mod chemistry;
mod common;
mod comparison;
mod crystal;
mod data_commands;
mod interactions;
mod mapped_compare;
mod sasa_tables;
mod sse;
mod structure;
mod trajectory;
mod validation;
mod validation_builtin;
