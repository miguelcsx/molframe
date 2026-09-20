//! Python's named facade over the typed native workflow.

use crate::bindings::{PyContactTable, PyStructure, coordinates};
use numpy::{IntoPyArray, PyArrayMethods, PyReadonlyArray2};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT_WORKFLOW_ID: AtomicU64 = AtomicU64::new(1);

#[derive(Clone, Copy, Debug)]
enum NodeValue {
    Coordinates(molframe::Node<Vec<[f32; 3]>>),
    Float(molframe::Node<f64>),
    Centroid(molframe::Node<Option<[f64; 3]>>),
    Matrix(molframe::Node<molframe::geometry::DistanceMatrix>),
    Structure(molframe::Node<molframe::Structure>),
    Contacts(molframe::Node<molframe::analysis::ContactTable>),
}

#[derive(Clone, Copy, Debug)]
#[pyclass(name = "WorkflowNode", frozen, skip_from_py_object)]
pub(crate) struct PyWorkflowNode {
    workflow_id: u64,
    value: NodeValue,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum InputKind {
    Coordinates,
    Float,
    Structure,
}

#[derive(Clone, Debug)]
enum OutputValue {
    Float(molframe::Output<f64>),
    Centroid(molframe::Output<Option<[f64; 3]>>),
    Matrix(molframe::Output<molframe::geometry::DistanceMatrix>),
    Structure(molframe::Output<molframe::Structure>),
    Contacts(molframe::Output<molframe::analysis::ContactTable>),
}

#[derive(Debug)]
#[pyclass(name = "Workflow", skip_from_py_object)]
pub(crate) struct PyWorkflow {
    id: u64,
    inner: molframe::Workflow,
    inputs: Vec<(Box<str>, InputKind)>,
    outputs: Vec<OutputValue>,
}

#[pymethods]
impl PyWorkflow {
    #[new]
    fn new() -> Self {
        Self {
            id: NEXT_WORKFLOW_ID.fetch_add(1, Ordering::Relaxed),
            inner: molframe::Workflow::new(),
            inputs: Vec::new(),
            outputs: Vec::new(),
        }
    }

    #[pyo3(signature = (name, *, kind, cost="borrow"))]
    fn input(&mut self, name: &str, kind: &str, cost: &str) -> PyResult<PyWorkflowNode> {
        let cost = parse_cost(cost)?;
        let (value, input_kind) = match kind {
            "coordinates" => (
                NodeValue::Coordinates(
                    self.inner
                        .input::<Vec<[f32; 3]>>(name, cost)
                        .map_err(workflow_error)?,
                ),
                InputKind::Coordinates,
            ),
            "float" => (
                NodeValue::Float(
                    self.inner
                        .input::<f64>(name, cost)
                        .map_err(workflow_error)?,
                ),
                InputKind::Float,
            ),
            "structure" => (
                NodeValue::Structure(
                    self.inner
                        .input::<molframe::Structure>(name, cost)
                        .map_err(workflow_error)?,
                ),
                InputKind::Structure,
            ),
            _ => {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "kind must be 'coordinates', 'float', or 'structure'",
                ));
            }
        };
        self.inputs.push((name.into(), input_kind));
        Ok(PyWorkflowNode {
            workflow_id: self.id,
            value,
        })
    }

    fn centroid(&mut self, coordinates: &PyWorkflowNode) -> PyResult<PyWorkflowNode> {
        self.ensure_node(coordinates)?;
        let NodeValue::Coordinates(input) = coordinates.value else {
            return Err(node_type("coordinates"));
        };
        let node = self
            .inner
            .map(
                input,
                molframe::OperationMetadata::new("geometry.centroid", molframe::Cost::Borrow)
                    .with_cse_key("geometry.centroid"),
                |coordinates, _| Ok(molframe::geometry::centroid(coordinates)),
            )
            .map_err(workflow_error)?;
        Ok(PyWorkflowNode {
            workflow_id: self.id,
            value: NodeValue::Centroid(node),
        })
    }

    fn rmsd(
        &mut self,
        mobile: &PyWorkflowNode,
        reference: &PyWorkflowNode,
    ) -> PyResult<PyWorkflowNode> {
        self.ensure_node(mobile)?;
        self.ensure_node(reference)?;
        let (NodeValue::Coordinates(mobile), NodeValue::Coordinates(reference)) =
            (mobile.value, reference.value)
        else {
            return Err(node_type("coordinates"));
        };
        let node = self
            .inner
            .map2(
                mobile,
                reference,
                molframe::OperationMetadata::new("geometry.rmsd", molframe::Cost::Borrow)
                    .with_cse_key("geometry.rmsd"),
                |mobile, reference, _| {
                    molframe::geometry::rmsd(mobile, reference).map_err(|error| {
                        molframe::WorkflowError::Operation {
                            operation: "geometry.rmsd",
                            message: format!("{error:?}").into(),
                        }
                    })
                },
            )
            .map_err(workflow_error)?;
        Ok(PyWorkflowNode {
            workflow_id: self.id,
            value: NodeValue::Float(node),
        })
    }

    fn distance_matrix(&mut self, coordinates: &PyWorkflowNode) -> PyResult<PyWorkflowNode> {
        self.ensure_node(coordinates)?;
        let NodeValue::Coordinates(input) = coordinates.value else {
            return Err(node_type("coordinates"));
        };
        let node =
            self.inner
                .map(
                    input,
                    molframe::OperationMetadata::new(
                        "geometry.distance_matrix",
                        molframe::Cost::Materialize,
                    )
                    .with_cse_key("geometry.distance_matrix"),
                    |coordinates, context| {
                        molframe::geometry::distance_matrix_with_context(coordinates, context)
                            .map_err(|error| molframe::WorkflowError::Operation {
                                operation: "geometry.distance_matrix",
                                message: error.to_string().into(),
                            })
                    },
                )
                .map_err(workflow_error)?;
        Ok(PyWorkflowNode {
            workflow_id: self.id,
            value: NodeValue::Matrix(node),
        })
    }

    #[pyo3(signature = (structure, cutoff))]
    fn atom_contacts(
        &mut self,
        structure: &PyWorkflowNode,
        cutoff: f32,
    ) -> PyResult<PyWorkflowNode> {
        self.ensure_node(structure)?;
        let NodeValue::Structure(input) = structure.value else {
            return Err(node_type("structure"));
        };
        let cse_key = format!(
            "analysis.atom_contacts:{cutoff:08x}",
            cutoff = cutoff.to_bits()
        );
        let node = self
            .inner
            .map(
                input,
                molframe::OperationMetadata::new(
                    "analysis.atom_contacts",
                    molframe::Cost::Materialize,
                )
                .with_cse_key(cse_key)
                .with_spatial_index(),
                move |structure, context| {
                    molframe::analysis::atom_contacts(
                        structure.engine(),
                        cutoff,
                        molframe::spatial::SpatialBackend::Auto,
                        context,
                    )
                    .map_err(|error| molframe::WorkflowError::Operation {
                        operation: "analysis.atom_contacts",
                        message: error.to_string().into(),
                    })
                },
            )
            .map_err(workflow_error)?;
        Ok(PyWorkflowNode {
            workflow_id: self.id,
            value: NodeValue::Contacts(node),
        })
    }

    fn output(&mut self, name: &str, node: &PyWorkflowNode) -> PyResult<()> {
        self.ensure_node(node)?;
        let output = match node.value {
            NodeValue::Float(node) => {
                OutputValue::Float(self.inner.output(name, node).map_err(workflow_error)?)
            }
            NodeValue::Centroid(node) => {
                OutputValue::Centroid(self.inner.output(name, node).map_err(workflow_error)?)
            }
            NodeValue::Matrix(node) => {
                OutputValue::Matrix(self.inner.output(name, node).map_err(workflow_error)?)
            }
            NodeValue::Structure(node) => {
                OutputValue::Structure(self.inner.output(name, node).map_err(workflow_error)?)
            }
            NodeValue::Contacts(node) => {
                OutputValue::Contacts(self.inner.output(name, node).map_err(workflow_error)?)
            }
            NodeValue::Coordinates(_) => {
                return Err(pyo3::exceptions::PyTypeError::new_err(
                    "raw coordinate inputs are not publishable outputs",
                ));
            }
        };
        self.outputs.push(output);
        Ok(())
    }

    fn compile(&self) -> PyResult<PyCompiledWorkflow> {
        Ok(PyCompiledWorkflow {
            inner: self.inner.compile().map_err(workflow_error)?,
            inputs: self.inputs.clone(),
            outputs: self.outputs.clone(),
        })
    }
}

impl PyWorkflow {
    fn ensure_node(&self, node: &PyWorkflowNode) -> PyResult<()> {
        if node.workflow_id == self.id {
            Ok(())
        } else {
            Err(pyo3::exceptions::PyValueError::new_err(
                "workflow node belongs to another Workflow",
            ))
        }
    }
}

#[derive(Clone, Debug)]
#[pyclass(name = "CompiledWorkflow", frozen, skip_from_py_object)]
pub(crate) struct PyCompiledWorkflow {
    inner: molframe::CompiledWorkflow,
    inputs: Vec<(Box<str>, InputKind)>,
    outputs: Vec<OutputValue>,
}

#[pymethods]
impl PyCompiledWorkflow {
    #[pyo3(signature = (values, *, copy=false))]
    fn run<'py>(
        &self,
        py: Python<'py>,
        values: &Bound<'py, PyDict>,
        copy: bool,
    ) -> PyResult<Bound<'py, PyDict>> {
        let mut inputs = molframe::WorkflowInputs::new();
        for (name, kind) in &self.inputs {
            let Some(value) = values.get_item(name.as_ref())? else {
                return Err(pyo3::exceptions::PyKeyError::new_err(name.to_string()));
            };
            match kind {
                InputKind::Coordinates => {
                    if !copy {
                        return Err(pyo3::exceptions::PyValueError::new_err(
                            "persistent coordinate inputs require copy=True",
                        ));
                    }
                    let array = value.extract::<PyReadonlyArray2<'_, f32>>()?;
                    inputs.insert(name.clone(), coordinates(&array)?.to_vec());
                }
                InputKind::Float => inputs.insert(name.clone(), value.extract::<f64>()?),
                InputKind::Structure => {
                    let structure = value.extract::<PyRef<'_, PyStructure>>()?;
                    inputs.insert(name.clone(), structure.inner.clone());
                }
            }
        }
        let results = py
            .detach(|| {
                self.inner
                    .run(&inputs, &molframe::ExecutionContext::default())
            })
            .map_err(workflow_error)?;
        let output = PyDict::new(py);
        for handle in &self.outputs {
            match handle {
                OutputValue::Float(handle) => {
                    output.set_item(handle.name(), results.get(handle).map_err(workflow_error)?)?;
                }
                OutputValue::Centroid(handle) => {
                    output.set_item(
                        handle.name(),
                        results.get(handle).map_err(workflow_error)?.as_ref(),
                    )?;
                }
                OutputValue::Matrix(handle) => {
                    let matrix = results.get(handle).map_err(workflow_error)?;
                    let array = matrix.as_slice().to_vec().into_pyarray(py);
                    output.set_item(
                        handle.name(),
                        array.reshape((matrix.rows(), matrix.columns()))?,
                    )?;
                }
                OutputValue::Structure(handle) => {
                    output.set_item(
                        handle.name(),
                        PyStructure::new(results.get(handle).map_err(workflow_error)?.clone()),
                    )?;
                }
                OutputValue::Contacts(handle) => {
                    let table = results.get_shared(handle).map_err(workflow_error)?;
                    output.set_item(handle.name(), PyContactTable::from_shared(py, table)?)?;
                }
            }
        }
        Ok(output)
    }

    fn explain<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let explanation = self.inner.explain();
        let result = PyDict::new(py);
        result.set_item("logical_node_count", explanation.logical_node_count)?;
        result.set_item("physical_node_count", explanation.physical_nodes.len())?;
        result.set_item(
            "common_subexpressions_eliminated",
            explanation.common_subexpressions_eliminated,
        )?;
        result.set_item("dead_nodes_eliminated", explanation.dead_nodes_eliminated)?;
        result.set_item("estimated_peak_bytes", explanation.estimated_peak_bytes)?;
        result.set_item("spatial_consumers", explanation.spatial_consumers)?;
        let nodes = PyList::empty(py);
        for node in &explanation.physical_nodes {
            let item = PyDict::new(py);
            item.set_item("ordinal", node.ordinal)?;
            item.set_item("operation", node.operation)?;
            item.set_item("cost", node.cost.to_string())?;
            item.set_item("estimated_output_bytes", node.estimated_output_bytes)?;
            nodes.append(item)?;
        }
        result.set_item("nodes", nodes)?;
        Ok(result)
    }
}

fn parse_cost(value: &str) -> PyResult<molframe::Cost> {
    match value {
        "borrow" => Ok(molframe::Cost::Borrow),
        "adopt" => Ok(molframe::Cost::Adopt),
        "decode" => Ok(molframe::Cost::Decode),
        "copy" => Ok(molframe::Cost::Copy),
        "materialize" => Ok(molframe::Cost::Materialize),
        _ => Err(pyo3::exceptions::PyValueError::new_err(
            "cost must be borrow, adopt, decode, copy, or materialize",
        )),
    }
}

fn node_type(expected: &str) -> PyErr {
    pyo3::exceptions::PyTypeError::new_err(format!("workflow node must contain {expected}"))
}

fn workflow_error(error: impl std::fmt::Display) -> PyErr {
    pyo3::exceptions::PyValueError::new_err(error.to_string())
}
