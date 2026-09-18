//! Mechanical bindings for native out-of-core provider contracts.

mod chunks;
mod metadata;
mod registration;
mod sources;

pub(crate) use registration::register;
