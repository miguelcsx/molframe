//! Typed, textual and reusable structural selection plans.

#![forbid(unsafe_code)]

mod annotation;
mod ast;
mod builder;
mod connectivity;
mod eval;
mod expand;
mod glob;
mod lexer;
mod macros;
mod parser;
mod plan;
mod predicate;
mod predicate_pattern;
mod spatial;

pub use builder::{Builder, ColumnBuilder, col};
pub use eval::{Evaluation, Groups, Query};
pub use plan::{LogicalPlan, PhysicalQuery};
pub use spatial::{GeometricRequest, SpatialResolver};
