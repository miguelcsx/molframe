//! Graph construction and the public typed handles.

use molframe_core::ExecutionContext;
use std::any::{Any, TypeId, type_name};
use std::fmt;
use std::marker::PhantomData;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use thiserror::Error;

pub(crate) type Value = Arc<dyn Any + Send + Sync>;
pub(crate) type Kernel =
    Arc<dyn Fn(&[Value], &ExecutionContext) -> WorkflowResult<Value> + Send + Sync>;

static NEXT_GRAPH_ID: AtomicU64 = AtomicU64::new(1);

/// The data-movement cost attached to an input or operation.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Cost {
    /// Reuses storage owned elsewhere for the duration of execution.
    Borrow,
    /// Retains compatible owner-backed storage without copying it.
    Adopt,
    /// Converts an encoded representation into native storage.
    Decode,
    /// Duplicates an existing native representation.
    Copy,
    /// Gathers or computes a new representation on explicit request.
    Materialize,
}

impl fmt::Display for Cost {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Borrow => "borrow",
            Self::Adopt => "adopt",
            Self::Decode => "decode",
            Self::Copy => "copy",
            Self::Materialize => "materialize",
        })
    }
}

/// Static information used by compilation, explanation and memory admission.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OperationMetadata {
    /// Canonical operation name shown by `explain`.
    pub name: &'static str,
    /// Stable semantic key. Equal keys and equal dependencies are CSE candidates.
    pub cse_key: Option<Box<str>>,
    /// Dominant data-movement cost.
    pub cost: Cost,
    /// Conservative bytes retained by the operation result.
    pub estimated_output_bytes: usize,
    /// Whether the operation needs a reusable spatial index.
    pub uses_spatial_index: bool,
}

impl OperationMetadata {
    /// Metadata for an operation without a semantic CSE key.
    #[must_use]
    pub const fn new(name: &'static str, cost: Cost) -> Self {
        Self {
            name,
            cse_key: None,
            cost,
            estimated_output_bytes: 0,
            uses_spatial_index: false,
        }
    }

    /// Adds a stable key used for common-subexpression elimination.
    #[must_use]
    pub fn with_cse_key(mut self, key: impl Into<Box<str>>) -> Self {
        self.cse_key = Some(key.into());
        self
    }

    /// Adds a conservative retained-size estimate.
    #[must_use]
    pub const fn with_estimated_output_bytes(mut self, bytes: usize) -> Self {
        self.estimated_output_bytes = bytes;
        self
    }

    /// Marks an operation as a spatial-index consumer.
    #[must_use]
    pub const fn with_spatial_index(mut self) -> Self {
        self.uses_spatial_index = true;
        self
    }
}

/// A typed immutable expression node.
pub struct Node<T> {
    pub(crate) graph: u64,
    pub(crate) index: usize,
    pub(crate) marker: PhantomData<fn() -> T>,
}

impl<T> Copy for Node<T> {}

impl<T> Clone for Node<T> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<T> fmt::Debug for Node<T> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("Node")
            .field("graph", &self.graph)
            .field("index", &self.index)
            .field("type", &type_name::<T>())
            .finish()
    }
}

/// A typed named workflow input.
pub type Input<T> = Node<T>;

/// A typed named workflow output handle.
#[derive(Clone, Debug)]
pub struct Output<T> {
    pub(crate) node: Node<T>,
    pub(crate) name: Box<str>,
}

impl<T> Output<T> {
    /// Output name used in named result lookup.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }
}

#[derive(Clone)]
pub(crate) struct GraphNode {
    pub(crate) type_id: TypeId,
    pub(crate) type_name: &'static str,
    pub(crate) dependencies: Vec<usize>,
    pub(crate) expected_inputs: Vec<TypeId>,
    pub(crate) input_name: Option<Box<str>>,
    pub(crate) metadata: OperationMetadata,
    pub(crate) kernel: Option<Kernel>,
}

impl fmt::Debug for GraphNode {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("GraphNode")
            .field("type_name", &self.type_name)
            .field("dependencies", &self.dependencies)
            .field("input_name", &self.input_name)
            .field("metadata", &self.metadata)
            .finish_non_exhaustive()
    }
}

/// A symbolic computation graph under construction.
#[derive(Clone, Debug)]
pub struct Workflow {
    pub(crate) id: u64,
    pub(crate) nodes: Vec<GraphNode>,
    pub(crate) outputs: Vec<(Box<str>, usize, TypeId)>,
}

impl Default for Workflow {
    fn default() -> Self {
        Self::new()
    }
}

impl Workflow {
    /// Starts an empty workflow.
    #[must_use]
    pub fn new() -> Self {
        Self {
            id: NEXT_GRAPH_ID.fetch_add(1, Ordering::Relaxed),
            nodes: Vec::new(),
            outputs: Vec::new(),
        }
    }

    /// Declares a typed named input.
    ///
    /// # Errors
    ///
    /// Returns an error when the name is empty or already used by an input.
    pub fn input<T: Send + Sync + 'static>(
        &mut self,
        name: impl Into<Box<str>>,
        cost: Cost,
    ) -> Result<Input<T>, WorkflowBuildError> {
        let name = name.into();
        if name.is_empty() {
            return Err(WorkflowBuildError::EmptyName);
        }
        if self
            .nodes
            .iter()
            .any(|node| node.input_name.as_deref() == Some(&name))
        {
            return Err(WorkflowBuildError::DuplicateInput(name));
        }
        let index = self.nodes.len();
        self.nodes.push(GraphNode {
            type_id: TypeId::of::<T>(),
            type_name: type_name::<T>(),
            dependencies: Vec::new(),
            expected_inputs: Vec::new(),
            input_name: Some(name),
            metadata: OperationMetadata::new("input", cost),
            kernel: None,
        });
        Ok(self.node(index))
    }

    /// Adds one typed unary operation.
    ///
    /// # Errors
    ///
    /// Returns an error when `input` belongs to another workflow.
    pub fn map<A, T, F>(
        &mut self,
        input: Node<A>,
        metadata: OperationMetadata,
        kernel: F,
    ) -> Result<Node<T>, WorkflowBuildError>
    where
        A: Send + Sync + 'static,
        T: Send + Sync + 'static,
        F: Fn(&A, &ExecutionContext) -> WorkflowResult<T> + Send + Sync + 'static,
    {
        self.ensure_graph(input.graph)?;
        let erased: Kernel = Arc::new(move |values, context| {
            let value = downcast::<A>(&values[0])?;
            kernel(value, context).map(|result| Arc::new(result) as Value)
        });
        Ok(self.push_operation::<T>(vec![input.index], vec![TypeId::of::<A>()], metadata, erased))
    }

    /// Adds one typed binary operation.
    ///
    /// # Errors
    ///
    /// Returns an error when either input belongs to another workflow.
    pub fn map2<A, B, T, F>(
        &mut self,
        first: Node<A>,
        second: Node<B>,
        metadata: OperationMetadata,
        kernel: F,
    ) -> Result<Node<T>, WorkflowBuildError>
    where
        A: Send + Sync + 'static,
        B: Send + Sync + 'static,
        T: Send + Sync + 'static,
        F: Fn(&A, &B, &ExecutionContext) -> WorkflowResult<T> + Send + Sync + 'static,
    {
        self.ensure_graph(first.graph)?;
        self.ensure_graph(second.graph)?;
        let erased: Kernel = Arc::new(move |values, context| {
            let first = downcast::<A>(&values[0])?;
            let second = downcast::<B>(&values[1])?;
            kernel(first, second, context).map(|result| Arc::new(result) as Value)
        });
        Ok(self.push_operation::<T>(
            vec![first.index, second.index],
            vec![TypeId::of::<A>(), TypeId::of::<B>()],
            metadata,
            erased,
        ))
    }

    /// Publishes a typed node under a stable output name.
    ///
    /// # Errors
    ///
    /// Returns an error for an empty/duplicate name or a foreign node.
    pub fn output<T: Send + Sync + 'static>(
        &mut self,
        name: impl Into<Box<str>>,
        node: Node<T>,
    ) -> Result<Output<T>, WorkflowBuildError> {
        self.ensure_graph(node.graph)?;
        let name = name.into();
        if name.is_empty() {
            return Err(WorkflowBuildError::EmptyName);
        }
        if self
            .outputs
            .iter()
            .any(|(existing, _, _)| existing == &name)
        {
            return Err(WorkflowBuildError::DuplicateOutput(name));
        }
        self.outputs
            .push((name.clone(), node.index, TypeId::of::<T>()));
        Ok(Output { node, name })
    }

    fn node<T>(&self, index: usize) -> Node<T> {
        Node {
            graph: self.id,
            index,
            marker: PhantomData,
        }
    }

    fn push_operation<T: Send + Sync + 'static>(
        &mut self,
        dependencies: Vec<usize>,
        expected_inputs: Vec<TypeId>,
        metadata: OperationMetadata,
        kernel: Kernel,
    ) -> Node<T> {
        let index = self.nodes.len();
        self.nodes.push(GraphNode {
            type_id: TypeId::of::<T>(),
            type_name: type_name::<T>(),
            dependencies,
            expected_inputs,
            input_name: None,
            metadata,
            kernel: Some(kernel),
        });
        self.node(index)
    }

    fn ensure_graph(&self, graph: u64) -> Result<(), WorkflowBuildError> {
        if graph == self.id {
            Ok(())
        } else {
            Err(WorkflowBuildError::ForeignNode)
        }
    }
}

fn downcast<T: Send + Sync + 'static>(value: &Value) -> WorkflowResult<&T> {
    value
        .downcast_ref::<T>()
        .ok_or(WorkflowError::InternalTypeMismatch {
            expected: type_name::<T>(),
        })
}

/// Graph-construction errors detected before compilation.
#[derive(Clone, Debug, Error, PartialEq, Eq)]
pub enum WorkflowBuildError {
    /// Names must be non-empty.
    #[error("workflow names must not be empty")]
    EmptyName,
    /// An input name was declared twice.
    #[error("duplicate workflow input: {0}")]
    DuplicateInput(Box<str>),
    /// An output name was declared twice.
    #[error("duplicate workflow output: {0}")]
    DuplicateOutput(Box<str>),
    /// A node from another workflow was supplied.
    #[error("node belongs to another workflow")]
    ForeignNode,
}

/// Workflow compilation or execution failure.
#[derive(Debug, Error)]
pub enum WorkflowError {
    /// The graph has no named outputs.
    #[error("workflow has no outputs")]
    NoOutputs,
    /// A dependency index is invalid.
    #[error("node {node} references missing dependency {dependency}")]
    MissingDependency {
        /// Node being validated.
        node: usize,
        /// Missing dependency index.
        dependency: usize,
    },
    /// A dependency type does not match the operation signature.
    #[error("node {node} has an incompatible dependency type")]
    TypeMismatch {
        /// Node being validated.
        node: usize,
    },
    /// A cycle prevents stable topological execution.
    #[error("workflow graph contains a cycle")]
    Cycle,
    /// A required named input was not supplied.
    #[error("missing workflow input: {0}")]
    MissingInput(Box<str>),
    /// A named input has the wrong concrete type.
    #[error("workflow input {name:?} has the wrong type; expected {expected}")]
    InputType {
        /// Input name.
        name: Box<str>,
        /// Expected Rust type.
        expected: &'static str,
    },
    /// Execution was cancelled between nodes.
    #[error("workflow execution cancelled")]
    Cancelled,
    /// Static output-lifetime planning exceeds the execution budget.
    #[error("workflow needs an estimated {estimated} bytes but the context permits {budget}")]
    MemoryRefused {
        /// Conservative workflow estimate.
        estimated: usize,
        /// Context memory ceiling.
        budget: usize,
    },
    /// A kernel rejected its input or could not produce its result.
    #[error("operation {operation:?} failed: {message}")]
    Operation {
        /// Operation name.
        operation: &'static str,
        /// Kernel error.
        message: Box<str>,
    },
    /// Internal type validation disagreed with the typed builder.
    #[error("internal workflow type mismatch; expected {expected}")]
    InternalTypeMismatch {
        /// Expected concrete Rust type.
        expected: &'static str,
    },
    /// A named result was requested with the wrong type.
    #[error("workflow result {name:?} has the wrong requested type")]
    ResultType {
        /// Result name.
        name: Box<str>,
    },
    /// A result handle does not belong to this execution.
    #[error("workflow result is unavailable: {0}")]
    MissingResult(Box<str>),
}

/// Result type used by workflow kernels.
pub type WorkflowResult<T> = Result<T, WorkflowError>;
