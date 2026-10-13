//! Detailed projections over native validation kernels.
//!
//! `validate --checks a,b,c` runs each named check against one structure and
//! prints one row per finding. A check that embodies a scientific definition
//! refuses to run without its explicit parameters.

mod basic;
mod pipeline;
mod reference;
mod structural;
mod types;

pub(crate) use pipeline::validate;
use types::{Row, Subject, ValidationOptions};
