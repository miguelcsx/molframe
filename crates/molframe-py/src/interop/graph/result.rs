//! One-call native graph materialisation and framework adapters.

use super::{PyEdgeFeature, PyGraphOptions, PyNodeFeature};
use crate::structure::PyStructure;
use molframe::Graph;
use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray2};
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

const NODE_DATA_KEY: &str = "x";
const EDGE_DATA_KEY: &str = "edge_attr";

#[pyclass(name = "Graph", frozen, skip_from_py_object)]
pub(crate) struct PyGraph {
    node_count: usize,
    edge_count: usize,
    edge_index: Py<PyArray2<i64>>,
    node_features: Py<PyArray2<f32>>,
    edge_features: Py<PyArray2<f32>>,
    node_feature_schema: Vec<PyNodeFeature>,
    edge_feature_schema: Vec<PyEdgeFeature>,
    cost: crate::extensions::PyExportCost,
}

#[pymethods]
impl PyStructure {
    fn graph(&self, py: Python<'_>, options: &PyGraphOptions) -> PyResult<PyGraph> {
        build_graph(py, self, options)
    }
}

#[pyfunction(name = "graph")]
pub(crate) fn build_graph(
    py: Python<'_>,
    structure: &PyStructure,
    options: &PyGraphOptions,
) -> PyResult<PyGraph> {
    molframe::graph(
        structure.structure(),
        &options.inner,
        &crate::core::execution::default_context(),
    )
    .map_err(|error| crate::errors::graph_error(&error))
    .and_then(|graph| PyGraph::new(py, graph))
}

#[pymethods]
impl PyGraph {
    #[getter]
    const fn node_count(&self) -> usize {
        self.node_count
    }

    #[getter]
    const fn edge_count(&self) -> usize {
        self.edge_count
    }

    #[getter]
    fn edge_index(&self, py: Python<'_>) -> Py<PyArray2<i64>> {
        self.edge_index.clone_ref(py)
    }

    #[getter]
    fn node_features(&self, py: Python<'_>) -> Py<PyArray2<f32>> {
        self.node_features.clone_ref(py)
    }

    #[getter]
    fn edge_features(&self, py: Python<'_>) -> Py<PyArray2<f32>> {
        self.edge_features.clone_ref(py)
    }

    #[getter]
    fn node_feature_schema(&self) -> Vec<PyNodeFeature> {
        self.node_feature_schema.clone()
    }

    #[getter]
    fn edge_feature_schema(&self) -> Vec<PyEdgeFeature> {
        self.edge_feature_schema.clone()
    }

    #[getter]
    fn cost(&self) -> crate::extensions::PyExportCost {
        self.cost
    }

    fn to_pyg(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let torch = py.import("torch")?;
        let data = py.import("torch_geometric.data")?.getattr("Data")?;
        let kwargs = PyDict::new(py);
        kwargs.set_item(
            NODE_DATA_KEY,
            torch.call_method1("from_numpy", (self.node_features.bind(py),))?,
        )?;
        kwargs.set_item(
            "edge_index",
            torch.call_method1("from_numpy", (self.edge_index.bind(py),))?,
        )?;
        kwargs.set_item(
            EDGE_DATA_KEY,
            torch.call_method1("from_numpy", (self.edge_features.bind(py),))?,
        )?;
        data.call((), Some(&kwargs)).map(Bound::unbind)
    }

    fn to_dgl(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let torch = py.import("torch")?;
        let dgl = py.import("dgl")?;
        let indices = self.edge_index.bind(py);
        let source = torch.call_method1("from_numpy", (indices.get_item(0)?,))?;
        let target = torch.call_method1("from_numpy", (indices.get_item(1)?,))?;
        let kwargs = PyDict::new(py);
        kwargs.set_item("num_nodes", self.node_count)?;
        let graph = dgl
            .getattr("graph")?
            .call(((source, target),), Some(&kwargs))?;
        graph.getattr("ndata")?.set_item(
            NODE_DATA_KEY,
            torch.call_method1("from_numpy", (self.node_features.bind(py),))?,
        )?;
        graph.getattr("edata")?.set_item(
            EDGE_DATA_KEY,
            torch.call_method1("from_numpy", (self.edge_features.bind(py),))?,
        )?;
        Ok(graph.unbind())
    }
}

impl PyGraph {
    fn new(py: Python<'_>, graph: Graph) -> PyResult<Self> {
        let edge_index = matrix(
            graph.edge_index.into_vec(),
            2,
            graph.edge_count,
            "edge index",
        )?
        .into_pyarray(py)
        .unbind();
        let node_feature_schema = graph
            .node_features
            .features
            .iter()
            .copied()
            .map(Into::into)
            .collect();
        let node_features = matrix(
            graph.node_features.values.into_vec(),
            graph.node_features.rows,
            graph.node_features.columns,
            "node features",
        )?
        .into_pyarray(py)
        .unbind();
        let edge_feature_schema = graph
            .edge_features
            .features
            .iter()
            .copied()
            .map(Into::into)
            .collect();
        let edge_features = matrix(
            graph.edge_features.values.into_vec(),
            graph.edge_features.rows,
            graph.edge_features.columns,
            "edge features",
        )?
        .into_pyarray(py)
        .unbind();
        Ok(Self {
            node_count: graph.node_count,
            edge_count: graph.edge_count,
            edge_index,
            node_features,
            edge_features,
            node_feature_schema,
            edge_feature_schema,
            cost: crate::extensions::PyExportCost::Copy,
        })
    }
}

fn matrix<T>(values: Vec<T>, rows: usize, columns: usize, name: &str) -> PyResult<Array2<T>> {
    Array2::from_shape_vec((rows, columns), values)
        .map_err(|error| PyRuntimeError::new_err(format!("invalid {name} matrix: {error}")))
}

#[cfg(test)]
#[path = "result_tests.rs"]
mod tests;
