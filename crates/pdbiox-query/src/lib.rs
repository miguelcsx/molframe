//! Typed, textual and reusable structural selection plans.

#![forbid(unsafe_code)]

mod api;
mod execution;
mod language;

pub(crate) use api::builder;
pub(crate) use execution::{annotation, connectivity, plan, spatial};
pub(crate) use language::{
    ast, expand, glob, lexer, macros, model_pattern, parser, predicate, predicate_pattern,
};

pub use api::{Builder, ColumnBuilder, col};
pub use execution::{
    Evaluation, GeometricRequest, Groups, LogicalPlan, PhysicalQuery, Query, SpatialResolver,
};
