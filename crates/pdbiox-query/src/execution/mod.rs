//! Query binding and execution.

pub(crate) mod annotation;
pub(crate) mod connectivity;
pub(crate) mod eval;
pub(crate) mod plan;
pub(crate) mod spatial;

pub use eval::{Evaluation, Groups, Query};
pub use plan::{LogicalPlan, PhysicalQuery};
pub use spatial::{GeometricRequest, SpatialResolver};
