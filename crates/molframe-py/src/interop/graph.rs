//! Graph tensors for graph neural network libraries.

use crate::bindings::PyStructure;
use molframe::interop::{
    EdgeDirection, EdgeFeature, EdgeKind, GraphOptions, MissingFeaturePolicy, NodeFeature,
    NodeLevel,
};
use numpy::{PyArray1, PyArray2, PyArrayMethods};
use pyo3::prelude::*;

/// Node and edge tensors in the layout `PyG` and `DGL` consume.
///
/// `edge_index` is `(2, edges)` `int64`, all sources then all targets;
/// `node_features` is `(nodes, features)` and `edge_features` `(edges,
/// features)`, both `float32` with their columns named in request order.
/// Every array is the caller's own copy.
#[derive(Debug)]
#[pyclass(
    name = "Graph",
    frozen,
    skip_from_py_object,
    module = "molframe.interop"
)]
pub(crate) struct PyGraph {
    node_count: usize,
    edge_count: usize,
    edge_index: Py<PyArray2<i64>>,
    node_features: Py<PyArray2<f32>>,
    edge_features: Py<PyArray2<f32>>,
    node_feature_names: Vec<String>,
    edge_feature_names: Vec<String>,
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

    /// The feature each column of `node_features` holds.
    #[getter]
    fn node_feature_names(&self) -> Vec<String> {
        self.node_feature_names.clone()
    }

    /// The feature each column of `edge_features` holds.
    #[getter]
    fn edge_feature_names(&self) -> Vec<String> {
        self.edge_feature_names.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "Graph(nodes={}, edges={}, node_features={}, edge_features={})",
            self.node_count,
            self.edge_count,
            self.node_feature_names.len(),
            self.edge_feature_names.len()
        )
    }
}

fn parse_all<T>(names: Vec<String>) -> PyResult<Vec<T>>
where
    T: std::str::FromStr<Err = molframe::PolicyParseError>,
{
    names
        .into_iter()
        .map(|name| name.parse().map_err(crate::error::kernel))
        .collect()
}

/// Builds a graph of the structure with every decision stated.
///
/// `nodes` is `"atoms"` or `"residues"`; `edges` is `"bonds"`, `"contacts"`
/// or `"radius"` (each with `cutoff` in ångström) or `"k_nearest"` (with
/// `neighbors`); `direction` is `"undirected"`, `"symmetric"` or `"directed"`.
/// `missing` is `None` to refuse a graph with an absent feature value, or the
/// finite value to fill it with. `periodic` measures through the unit cell.
#[pyfunction]
#[pyo3(signature = (
    structure,
    *,
    nodes,
    edges,
    direction,
    cutoff=None,
    neighbors=None,
    node_features=Vec::new(),
    edge_features=Vec::new(),
    missing=None,
    backend="auto",
    periodic=false,
    context=None,
))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn graph(
    py: Python<'_>,
    structure: &PyStructure,
    nodes: &str,
    edges: &str,
    direction: &str,
    cutoff: Option<f32>,
    neighbors: Option<usize>,
    node_features: Vec<String>,
    edge_features: Vec<String>,
    missing: Option<f32>,
    backend: &str,
    periodic: bool,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<PyGraph> {
    let node_features: Vec<NodeFeature> = parse_all(node_features)?;
    let edge_features: Vec<EdgeFeature> = parse_all(edge_features)?;
    let options = GraphOptions {
        nodes: nodes.parse::<NodeLevel>().map_err(crate::error::kernel)?,
        edges: EdgeKind::from_parts(edges, cutoff, neighbors).map_err(crate::error::kernel)?,
        direction: direction
            .parse::<EdgeDirection>()
            .map_err(crate::error::kernel)?,
        node_features,
        edge_features,
        missing: missing.map_or(MissingFeaturePolicy::Error, MissingFeaturePolicy::Fill),
        backend: crate::backend::parse(backend)?,
        periodic,
    };
    let source = structure.inner.clone();
    let built = crate::execution::run(py, context, |context| {
        molframe::interop::graph(source.engine(), &options, context)
    })?
    .map_err(crate::error::kernel)?;
    let node_names: Vec<String> = built
        .node_features
        .features
        .iter()
        .map(|feature| feature.name().replace('-', "_"))
        .collect();
    let edge_names: Vec<String> = built
        .edge_features
        .features
        .iter()
        .map(|feature| feature.name().replace('-', "_"))
        .collect();
    let matrix = |values: Box<[f32]>, rows: usize, columns: usize| {
        PyArray1::from_vec(py, values.into_vec()).reshape([rows, columns])
    };
    Ok(PyGraph {
        node_count: built.node_count,
        edge_count: built.edge_count,
        edge_index: PyArray1::from_vec(py, built.edge_index.into_vec())
            .reshape([2, built.edge_count])?
            .unbind(),
        node_features: matrix(
            built.node_features.values,
            built.node_features.rows,
            built.node_features.columns,
        )?
        .unbind(),
        edge_features: matrix(
            built.edge_features.values,
            built.edge_features.rows,
            built.edge_features.columns,
        )?
        .unbind(),
        node_feature_names: node_names,
        edge_feature_names: edge_names,
    })
}
