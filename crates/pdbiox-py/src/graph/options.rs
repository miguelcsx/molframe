//! Typed graph policies with no stringly-typed dispatch.

use pdbiox::{
    EdgeDirection, EdgeFeature, EdgeKind, GraphOptions, MissingFeaturePolicy, NodeFeature,
    NodeLevel, SpatialBackend,
};
use pyo3::prelude::*;

macro_rules! fieldless_enum {
    ($python:literal, $name:ident { $($variant:ident),+ $(,)? }) => {
        #[pyclass(name = $python, frozen, eq, eq_int, from_py_object)]
        #[derive(Clone, Copy, Debug, PartialEq, Eq)]
        pub(crate) enum $name { $($variant),+ }
    };
}

fieldless_enum!("NodeLevel", PyNodeLevel { Atoms, Residues });
fieldless_enum!(
    "EdgeDirection",
    PyEdgeDirection {
        Undirected,
        Symmetric,
        Directed
    }
);
fieldless_enum!(
    "NodeFeature",
    PyNodeFeature {
        Element,
        FormalCharge,
        PartialCharge,
        BFactor,
        Occupancy,
        AtomCount,
        PositionX,
        PositionY,
        PositionZ,
    }
);
fieldless_enum!(
    "EdgeFeature",
    PyEdgeFeature {
        Distance,
        BondOrder
    }
);
fieldless_enum!(
    "SpatialBackend",
    PySpatialBackend {
        BruteForce,
        CellList,
        KdTree,
        NeighborList,
        Auto,
    }
);

#[pyclass(name = "EdgeKind", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyEdgeKind(EdgeKind);

#[pymethods]
impl PyEdgeKind {
    #[staticmethod]
    fn bonds() -> Self {
        Self(EdgeKind::Bonds)
    }

    #[staticmethod]
    fn contacts(cutoff: f32) -> Self {
        Self(EdgeKind::Contacts { cutoff })
    }

    #[staticmethod]
    fn radius(cutoff: f32) -> Self {
        Self(EdgeKind::Radius { cutoff })
    }

    #[staticmethod]
    fn k_nearest(neighbors: usize) -> Self {
        Self(EdgeKind::KNearest { neighbors })
    }
}

#[pyclass(name = "MissingFeaturePolicy", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyMissingFeaturePolicy(MissingFeaturePolicy);

#[pymethods]
impl PyMissingFeaturePolicy {
    #[staticmethod]
    fn error() -> Self {
        Self(MissingFeaturePolicy::Error)
    }

    #[staticmethod]
    fn fill(value: f32) -> Self {
        Self(MissingFeaturePolicy::Fill(value))
    }
}

#[pyclass(name = "GraphOptions", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyGraphOptions {
    pub(super) inner: GraphOptions,
}

#[pymethods]
impl PyGraphOptions {
    #[new]
    #[pyo3(signature = (nodes, edges, *, direction=PyEdgeDirection::Symmetric, features=(Vec::new(), Vec::new()), missing=None, backend=PySpatialBackend::Auto, periodic=false))]
    fn new(
        nodes: PyNodeLevel,
        edges: PyEdgeKind,
        direction: PyEdgeDirection,
        features: (Vec<PyNodeFeature>, Vec<PyEdgeFeature>),
        missing: Option<PyMissingFeaturePolicy>,
        backend: PySpatialBackend,
        periodic: bool,
    ) -> Self {
        Self {
            inner: GraphOptions {
                nodes: nodes.into(),
                edges: edges.0,
                direction: direction.into(),
                node_features: features.0.into_iter().map(Into::into).collect(),
                edge_features: features.1.into_iter().map(Into::into).collect(),
                missing: missing.map_or(MissingFeaturePolicy::Error, |value| value.0),
                backend: backend.into(),
                periodic,
            },
        }
    }
}

impl From<PyNodeLevel> for NodeLevel {
    fn from(value: PyNodeLevel) -> Self {
        match value {
            PyNodeLevel::Atoms => Self::Atoms,
            PyNodeLevel::Residues => Self::Residues,
        }
    }
}

impl From<PyEdgeDirection> for EdgeDirection {
    fn from(value: PyEdgeDirection) -> Self {
        match value {
            PyEdgeDirection::Undirected => Self::Undirected,
            PyEdgeDirection::Symmetric => Self::Symmetric,
            PyEdgeDirection::Directed => Self::Directed,
        }
    }
}

impl From<PyNodeFeature> for NodeFeature {
    fn from(value: PyNodeFeature) -> Self {
        match value {
            PyNodeFeature::Element => Self::Element,
            PyNodeFeature::FormalCharge => Self::FormalCharge,
            PyNodeFeature::PartialCharge => Self::PartialCharge,
            PyNodeFeature::BFactor => Self::BFactor,
            PyNodeFeature::Occupancy => Self::Occupancy,
            PyNodeFeature::AtomCount => Self::AtomCount,
            PyNodeFeature::PositionX => Self::PositionX,
            PyNodeFeature::PositionY => Self::PositionY,
            PyNodeFeature::PositionZ => Self::PositionZ,
        }
    }
}

impl From<NodeFeature> for PyNodeFeature {
    fn from(value: NodeFeature) -> Self {
        match value {
            NodeFeature::Element => Self::Element,
            NodeFeature::FormalCharge => Self::FormalCharge,
            NodeFeature::PartialCharge => Self::PartialCharge,
            NodeFeature::BFactor => Self::BFactor,
            NodeFeature::Occupancy => Self::Occupancy,
            NodeFeature::AtomCount => Self::AtomCount,
            NodeFeature::PositionX => Self::PositionX,
            NodeFeature::PositionY => Self::PositionY,
            NodeFeature::PositionZ => Self::PositionZ,
        }
    }
}

impl From<PyEdgeFeature> for EdgeFeature {
    fn from(value: PyEdgeFeature) -> Self {
        match value {
            PyEdgeFeature::Distance => Self::Distance,
            PyEdgeFeature::BondOrder => Self::BondOrder,
        }
    }
}

impl From<EdgeFeature> for PyEdgeFeature {
    fn from(value: EdgeFeature) -> Self {
        match value {
            EdgeFeature::Distance => Self::Distance,
            EdgeFeature::BondOrder => Self::BondOrder,
        }
    }
}

impl From<PySpatialBackend> for SpatialBackend {
    fn from(value: PySpatialBackend) -> Self {
        match value {
            PySpatialBackend::BruteForce => Self::BruteForce,
            PySpatialBackend::CellList => Self::CellList,
            PySpatialBackend::KdTree => Self::KdTree,
            PySpatialBackend::NeighborList => Self::NeighborList,
            PySpatialBackend::Auto => Self::Auto,
        }
    }
}

impl From<SpatialBackend> for PySpatialBackend {
    fn from(value: SpatialBackend) -> Self {
        match value {
            SpatialBackend::BruteForce => Self::BruteForce,
            SpatialBackend::CellList => Self::CellList,
            SpatialBackend::KdTree => Self::KdTree,
            SpatialBackend::NeighborList => Self::NeighborList,
            _ => Self::Auto,
        }
    }
}
