//! Typed, reusable execution graphs for `MolFrame` kernels.
//!
//! The engine owns graph mechanics only. Domain crates and the public facade
//! add operations by closing over their statically dispatched kernels, keeping
//! this crate independent of `molframe` and preventing a dependency cycle.

#![forbid(unsafe_code)]

mod compile;
mod graph;
mod runtime;

pub use compile::{CompiledWorkflow, Explanation, PhysicalNode};
pub use graph::{
    Cost, Input, Node, OperationMetadata, Output, Workflow, WorkflowBuildError, WorkflowError,
    WorkflowResult,
};
pub use runtime::{WorkflowInputs, WorkflowResults};

#[cfg(test)]
#[path = "workflow_tests.rs"]
mod tests;
