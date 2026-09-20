//! Typed input ownership and deterministic execution.

use crate::compile::CompiledWorkflow;
use crate::graph::{Output, Value, WorkflowError};
use molframe_core::ExecutionContext;
use std::any::{TypeId, type_name};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

#[derive(Clone)]
struct InputValue {
    value: Value,
    type_id: TypeId,
}

/// Owned values supplied to a compiled workflow by input name.
#[derive(Clone, Default)]
pub struct WorkflowInputs {
    values: BTreeMap<Box<str>, InputValue>,
}

impl std::fmt::Debug for WorkflowInputs {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("WorkflowInputs")
            .field("names", &self.values.keys().collect::<Vec<_>>())
            .finish()
    }
}

impl WorkflowInputs {
    /// Starts an empty named input set.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            values: BTreeMap::new(),
        }
    }

    /// Inserts or replaces one typed value.
    pub fn insert<T: Send + Sync + 'static>(&mut self, name: impl Into<Box<str>>, value: T) {
        self.values.insert(
            name.into(),
            InputValue {
                value: Arc::new(value),
                type_id: TypeId::of::<T>(),
            },
        );
    }
}

/// Named, owner-backed values returned by one workflow execution.
#[derive(Clone, Default)]
pub struct WorkflowResults {
    graph_id: u64,
    values: BTreeMap<Box<str>, Value>,
}

impl std::fmt::Debug for WorkflowResults {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("WorkflowResults")
            .field("graph_id", &self.graph_id)
            .field("names", &self.values.keys().collect::<Vec<_>>())
            .finish()
    }
}

impl WorkflowResults {
    /// Returns a result through its typed output handle.
    ///
    /// # Errors
    ///
    /// Returns an error when the handle belongs to another graph, the result
    /// is missing or its requested type is inconsistent.
    pub fn get<T: Send + Sync + 'static>(&self, output: &Output<T>) -> Result<&T, WorkflowError> {
        if output.node.graph != self.graph_id {
            return Err(WorkflowError::MissingResult(output.name.clone()));
        }
        self.get_named(&output.name)
    }

    /// Clones the result's shared owner without copying its payload.
    ///
    /// This is the lifetime-safe path for FFI views that must outlive the
    /// [`WorkflowResults`] map.
    ///
    /// # Errors
    ///
    /// Returns the same handle, name, and type errors as [`Self::get`].
    pub fn get_shared<T: Send + Sync + 'static>(
        &self,
        output: &Output<T>,
    ) -> Result<Arc<T>, WorkflowError> {
        if output.node.graph != self.graph_id {
            return Err(WorkflowError::MissingResult(output.name.clone()));
        }
        let Some(value) = self.values.get(&output.name) else {
            return Err(WorkflowError::MissingResult(output.name.clone()));
        };
        Arc::clone(value)
            .downcast::<T>()
            .map_err(|_| WorkflowError::ResultType {
                name: output.name.clone(),
            })
    }

    /// Returns a result by its stable name and expected Rust type.
    ///
    /// # Errors
    ///
    /// Returns an error when the name is absent or the requested type differs.
    pub fn get_named<T: Send + Sync + 'static>(&self, name: &str) -> Result<&T, WorkflowError> {
        let Some(value) = self.values.get(name) else {
            return Err(WorkflowError::MissingResult(name.into()));
        };
        value
            .downcast_ref::<T>()
            .ok_or_else(|| WorkflowError::ResultType { name: name.into() })
    }

    /// Stable output names in lexical order.
    pub fn names(&self) -> impl ExactSizeIterator<Item = &str> {
        self.values.keys().map(AsRef::as_ref)
    }
}

pub(crate) fn execute(
    workflow: &CompiledWorkflow,
    inputs: &WorkflowInputs,
    context: &ExecutionContext,
) -> Result<WorkflowResults, WorkflowError> {
    if workflow.explanation.estimated_peak_bytes > context.memory_budget().bytes() {
        return Err(WorkflowError::MemoryRefused {
            estimated: workflow.explanation.estimated_peak_bytes,
            budget: context.memory_budget().bytes(),
        });
    }
    let _reservation = context
        .try_reserve(workflow.explanation.estimated_peak_bytes)
        .map_err(|_| WorkflowError::MemoryRefused {
            estimated: workflow.explanation.estimated_peak_bytes,
            budget: context.memory_budget().bytes(),
        })?;
    let output_nodes = workflow
        .outputs
        .iter()
        .map(|(_, node, _)| *node)
        .collect::<BTreeSet<_>>();
    let capacity = workflow
        .nodes
        .iter()
        .map(|node| node.index)
        .max()
        .map_or(0, |maximum| maximum + 1);
    let mut values = vec![None; capacity];
    let mut uses = workflow.use_counts.clone();
    for node in &workflow.nodes {
        if context.cancellation().is_cancelled() {
            return Err(WorkflowError::Cancelled);
        }
        let value = if let Some(name) = &node.input_name {
            let Some(input) = inputs.values.get(name) else {
                return Err(WorkflowError::MissingInput(name.clone()));
            };
            if input.type_id != node.type_id {
                return Err(WorkflowError::InputType {
                    name: name.clone(),
                    expected: node.type_name,
                });
            }
            Arc::clone(&input.value)
        } else {
            let arguments = node
                .dependencies
                .iter()
                .map(|dependency| {
                    values
                        .get(*dependency)
                        .and_then(Option::as_ref)
                        .cloned()
                        .ok_or(WorkflowError::InternalTypeMismatch {
                            expected: type_name::<()>(),
                        })
                })
                .collect::<Result<Vec<_>, _>>()?;
            let Some(kernel) = &node.kernel else {
                return Err(WorkflowError::InternalTypeMismatch {
                    expected: node.type_name,
                });
            };
            match kernel(&arguments, context) {
                Ok(value) => value,
                Err(error @ WorkflowError::Operation { .. }) => return Err(error),
                Err(error) => {
                    return Err(WorkflowError::Operation {
                        operation: node.operation,
                        message: error.to_string().into(),
                    });
                }
            }
        };
        values[node.index] = Some(value);
        for dependency in &node.dependencies {
            uses[*dependency] = uses[*dependency].saturating_sub(1);
            if uses[*dependency] == 0 && !output_nodes.contains(dependency) {
                values[*dependency] = None;
            }
        }
    }
    let mut results = BTreeMap::new();
    for (name, node, _) in &workflow.outputs {
        let Some(value) = values.get(*node).and_then(Option::as_ref) else {
            return Err(WorkflowError::MissingResult(name.clone()));
        };
        results.insert(name.clone(), Arc::clone(value));
    }
    Ok(WorkflowResults {
        graph_id: workflow.graph_id,
        values: results,
    })
}
