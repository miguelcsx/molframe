//! Validation, graph rewrites and the reusable physical plan.

use crate::graph::{GraphNode, Kernel, Workflow, WorkflowError};
use crate::runtime::{WorkflowInputs, WorkflowResults};
use molframe_core::ExecutionContext;
use std::any::TypeId;
use std::collections::{BTreeSet, HashMap, VecDeque};

/// One node in a compiled physical explanation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PhysicalNode {
    /// Stable physical execution ordinal.
    pub ordinal: usize,
    /// Original graph node index.
    pub node: usize,
    /// Canonical operation name.
    pub operation: &'static str,
    /// Canonical dependency node indices.
    pub dependencies: Vec<usize>,
    /// Dominant data-movement cost.
    pub cost: crate::Cost,
    /// Conservative retained result bytes.
    pub estimated_output_bytes: usize,
}

/// Logical and physical views of a compiled workflow.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Explanation {
    /// Nodes declared in the logical graph, before dead-node removal.
    pub logical_node_count: usize,
    /// Nodes retained by the physical graph.
    pub physical_nodes: Vec<PhysicalNode>,
    /// Logical nodes removed as common subexpressions.
    pub common_subexpressions_eliminated: usize,
    /// Logical nodes removed because no output reaches them.
    pub dead_nodes_eliminated: usize,
    /// Conservative peak retained bytes from output-lifetime planning.
    pub estimated_peak_bytes: usize,
    /// Number of retained operations that can share spatial planning state.
    pub spatial_consumers: usize,
}

#[derive(Clone)]
pub(crate) struct CompiledNode {
    pub(crate) index: usize,
    pub(crate) dependencies: Vec<usize>,
    pub(crate) input_name: Option<Box<str>>,
    pub(crate) type_id: TypeId,
    pub(crate) type_name: &'static str,
    pub(crate) operation: &'static str,
    pub(crate) kernel: Option<Kernel>,
}

impl std::fmt::Debug for CompiledNode {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("CompiledNode")
            .field("index", &self.index)
            .field("dependencies", &self.dependencies)
            .field("input_name", &self.input_name)
            .field("type_name", &self.type_name)
            .field("operation", &self.operation)
            .finish_non_exhaustive()
    }
}

/// A validated, optimized workflow reusable across input sets.
#[derive(Clone, Debug)]
pub struct CompiledWorkflow {
    pub(crate) graph_id: u64,
    pub(crate) nodes: Vec<CompiledNode>,
    pub(crate) outputs: Vec<(Box<str>, usize, TypeId)>,
    pub(crate) use_counts: Vec<usize>,
    pub(crate) explanation: Explanation,
}

impl Workflow {
    /// Validates and optimizes this graph into a reusable physical plan.
    ///
    /// Compilation performs type and cycle validation, stable topological
    /// ordering, common-subexpression elimination, dead-node elimination,
    /// output-lifetime planning and a conservative memory estimate.
    ///
    /// # Errors
    ///
    /// Returns an error for a graph with no outputs, an invalid dependency,
    /// an incompatible dependency type or a cycle.
    pub fn compile(&self) -> Result<CompiledWorkflow, WorkflowError> {
        if self.outputs.is_empty() {
            return Err(WorkflowError::NoOutputs);
        }
        validate_types(&self.nodes)?;
        let logical_order = stable_topological_order(&self.nodes)?;
        let aliases = common_subexpressions(&self.nodes, &logical_order);
        let canonical_outputs = self
            .outputs
            .iter()
            .map(|(name, node, type_id)| (name.clone(), aliases[*node], *type_id))
            .collect::<Vec<_>>();
        let active = reachable_nodes(&self.nodes, &aliases, &canonical_outputs);
        let physical_order = logical_order
            .iter()
            .copied()
            .filter(|node| active.contains(&aliases[*node]) && aliases[*node] == *node)
            .collect::<Vec<_>>();
        let nodes = physical_order
            .iter()
            .map(|index| {
                let node = &self.nodes[*index];
                CompiledNode {
                    index: *index,
                    dependencies: node
                        .dependencies
                        .iter()
                        .map(|dependency| aliases[*dependency])
                        .collect(),
                    input_name: node.input_name.clone(),
                    type_id: node.type_id,
                    type_name: node.type_name,
                    operation: node.metadata.name,
                    kernel: node.kernel.clone(),
                }
            })
            .collect::<Vec<_>>();
        let use_counts = use_counts(self.nodes.len(), &nodes, &canonical_outputs);
        let estimated_peak_bytes =
            estimate_peak(&self.nodes, &nodes, &canonical_outputs, &use_counts);
        let physical_nodes = nodes
            .iter()
            .enumerate()
            .map(|(ordinal, node)| PhysicalNode {
                ordinal,
                node: node.index,
                operation: node.operation,
                dependencies: node.dependencies.clone(),
                cost: self.nodes[node.index].metadata.cost,
                estimated_output_bytes: self.nodes[node.index].metadata.estimated_output_bytes,
            })
            .collect();
        let common_subexpressions_eliminated = aliases
            .iter()
            .enumerate()
            .filter(|(node, canonical)| *node != **canonical)
            .count();
        let dead_nodes_eliminated = self
            .nodes
            .len()
            .saturating_sub(nodes.len() + common_subexpressions_eliminated);
        let spatial_consumers = nodes
            .iter()
            .filter(|node| self.nodes[node.index].metadata.uses_spatial_index)
            .count();
        Ok(CompiledWorkflow {
            graph_id: self.id,
            nodes,
            outputs: canonical_outputs,
            use_counts,
            explanation: Explanation {
                logical_node_count: self.nodes.len(),
                physical_nodes,
                common_subexpressions_eliminated,
                dead_nodes_eliminated,
                estimated_peak_bytes,
                spatial_consumers,
            },
        })
    }
}

impl CompiledWorkflow {
    /// Returns both the logical summary and optimized physical plan.
    #[must_use]
    pub const fn explain(&self) -> &Explanation {
        &self.explanation
    }

    /// Executes the physical graph in stable topological order.
    ///
    /// # Errors
    ///
    /// Returns a missing/type input error, cancellation, memory refusal or the
    /// first kernel failure. No partial result is published.
    pub fn run(
        &self,
        inputs: &WorkflowInputs,
        context: &ExecutionContext,
    ) -> Result<WorkflowResults, WorkflowError> {
        crate::runtime::execute(self, inputs, context)
    }
}

fn validate_types(nodes: &[GraphNode]) -> Result<(), WorkflowError> {
    for (index, node) in nodes.iter().enumerate() {
        if node.dependencies.len() != node.expected_inputs.len() {
            return Err(WorkflowError::TypeMismatch { node: index });
        }
        for (dependency, expected) in node.dependencies.iter().zip(&node.expected_inputs) {
            let Some(actual) = nodes.get(*dependency) else {
                return Err(WorkflowError::MissingDependency {
                    node: index,
                    dependency: *dependency,
                });
            };
            if actual.type_id != *expected {
                return Err(WorkflowError::TypeMismatch { node: index });
            }
        }
    }
    Ok(())
}

fn stable_topological_order(nodes: &[GraphNode]) -> Result<Vec<usize>, WorkflowError> {
    let mut incoming = nodes
        .iter()
        .map(|node| node.dependencies.len())
        .collect::<Vec<_>>();
    let mut outgoing = vec![Vec::new(); nodes.len()];
    for (node, record) in nodes.iter().enumerate() {
        for dependency in &record.dependencies {
            outgoing[*dependency].push(node);
        }
    }
    let mut ready = incoming
        .iter()
        .enumerate()
        .filter_map(|(index, count)| (*count == 0).then_some(index))
        .collect::<BTreeSet<_>>();
    let mut order = Vec::with_capacity(nodes.len());
    while let Some(index) = ready.pop_first() {
        order.push(index);
        for dependent in &outgoing[index] {
            incoming[*dependent] -= 1;
            if incoming[*dependent] == 0 {
                ready.insert(*dependent);
            }
        }
    }
    if order.len() == nodes.len() {
        Ok(order)
    } else {
        Err(WorkflowError::Cycle)
    }
}

fn common_subexpressions(nodes: &[GraphNode], order: &[usize]) -> Vec<usize> {
    let mut aliases = (0..nodes.len()).collect::<Vec<_>>();
    let mut seen: HashMap<(Box<str>, Vec<usize>, TypeId), usize> = HashMap::new();
    for index in order {
        let node = &nodes[*index];
        let Some(key) = &node.metadata.cse_key else {
            continue;
        };
        let dependencies = node
            .dependencies
            .iter()
            .map(|dependency| aliases[*dependency])
            .collect::<Vec<_>>();
        let signature = (key.clone(), dependencies, node.type_id);
        if let Some(canonical) = seen.get(&signature) {
            aliases[*index] = *canonical;
        } else {
            seen.insert(signature, *index);
        }
    }
    aliases
}

fn reachable_nodes(
    nodes: &[GraphNode],
    aliases: &[usize],
    outputs: &[(Box<str>, usize, TypeId)],
) -> BTreeSet<usize> {
    let mut active = BTreeSet::new();
    let mut queue = outputs
        .iter()
        .map(|(_, node, _)| *node)
        .collect::<VecDeque<_>>();
    while let Some(index) = queue.pop_front() {
        if !active.insert(index) {
            continue;
        }
        for dependency in &nodes[index].dependencies {
            queue.push_back(aliases[*dependency]);
        }
    }
    active
}

fn use_counts(
    node_count: usize,
    nodes: &[CompiledNode],
    outputs: &[(Box<str>, usize, TypeId)],
) -> Vec<usize> {
    let mut counts = vec![0; node_count];
    for node in nodes {
        for dependency in &node.dependencies {
            counts[*dependency] += 1;
        }
    }
    for (_, node, _) in outputs {
        counts[*node] += 1;
    }
    counts
}

fn estimate_peak(
    logical: &[GraphNode],
    nodes: &[CompiledNode],
    outputs: &[(Box<str>, usize, TypeId)],
    initial_counts: &[usize],
) -> usize {
    let output_nodes = outputs
        .iter()
        .map(|(_, node, _)| *node)
        .collect::<BTreeSet<_>>();
    let mut counts = initial_counts.to_vec();
    let mut live = 0usize;
    let mut peak = 0usize;
    for node in nodes {
        live = live.saturating_add(logical[node.index].metadata.estimated_output_bytes);
        peak = peak.max(live);
        for dependency in &node.dependencies {
            counts[*dependency] = counts[*dependency].saturating_sub(1);
            if counts[*dependency] == 0 && !output_nodes.contains(dependency) {
                live = live.saturating_sub(logical[*dependency].metadata.estimated_output_bytes);
            }
        }
    }
    peak
}
