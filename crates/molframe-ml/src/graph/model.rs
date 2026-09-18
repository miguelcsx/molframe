//! Public graph schema and explicit construction policies.

use molframe_spatial::{SpatialBackend, SpatialError};

/// Rows represented as graph nodes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NodeLevel {
    /// One node per atom.
    Atoms,
    /// One node per residue, positioned at its finite-coordinate centroid.
    Residues,
}

/// Edge construction algorithm.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum EdgeKind {
    /// Explicit chemical bonds. Only atom graphs support this mode.
    Bonds,
    /// Fixed-radius non-bonded contacts.
    Contacts {
        /// Maximum node separation in angstroms.
        cutoff: f32,
    },
    /// Fixed-radius geometric graph including bonded pairs.
    Radius {
        /// Maximum node separation in angstroms.
        cutoff: f32,
    },
    /// Nearest neighbours per node.
    KNearest {
        /// Requested neighbours per node.
        neighbors: usize,
    },
}

/// Edge-index orientation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EdgeDirection {
    /// Store every relation once with its lower endpoint first.
    Undirected,
    /// Store both orientations of every relation.
    Symmetric,
    /// Preserve the query-to-neighbour direction for nearest-neighbour graphs.
    Directed,
}

/// Numeric node feature requested by name and type.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum NodeFeature {
    /// Atomic number; atom graphs only.
    Element,
    /// Formal charge assigned by structure chemistry.
    FormalCharge,
    /// Partial charge annotation.
    PartialCharge,
    /// Deposited temperature factor, or residue mean.
    BFactor,
    /// Deposited occupancy, or residue mean.
    Occupancy,
    /// Number of atoms represented by a node; residue graphs only.
    AtomCount,
    /// Cartesian x coordinate of the node position.
    PositionX,
    /// Cartesian y coordinate of the node position.
    PositionY,
    /// Cartesian z coordinate of the node position.
    PositionZ,
}

/// Numeric edge feature requested by name and type.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum EdgeFeature {
    /// Euclidean or minimum-image endpoint separation in angstroms.
    Distance,
    /// Numeric chemical bond order; explicit bond graphs only.
    BondOrder,
}

/// Required treatment when a requested value is absent.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MissingFeaturePolicy {
    /// Refuse the graph rather than fabricate a value.
    Error,
    /// Fill absent cells with this explicit finite value.
    Fill(f32),
}

/// Complete, explicit graph construction request.
#[derive(Clone, Debug, PartialEq)]
pub struct GraphOptions {
    /// Node granularity.
    pub nodes: NodeLevel,
    /// Edge algorithm and its parameters.
    pub edges: EdgeKind,
    /// Edge-index orientation.
    pub direction: EdgeDirection,
    /// Ordered node-feature columns.
    pub node_features: Vec<NodeFeature>,
    /// Ordered edge-feature columns.
    pub edge_features: Vec<EdgeFeature>,
    /// Missing-value behavior.
    pub missing: MissingFeaturePolicy,
    /// Spatial implementation used by geometric edges.
    pub backend: SpatialBackend,
    /// Whether geometric distances use the structure unit cell.
    pub periodic: bool,
}

/// Dense row-major feature matrix with typed ordered columns.
#[derive(Clone, Debug, PartialEq)]
pub struct FeatureMatrix<F> {
    /// Row-major `float32` values.
    pub values: Box<[f32]>,
    /// Matrix row count.
    pub rows: usize,
    /// Matrix column count.
    pub columns: usize,
    /// Feature descriptor for every column.
    pub features: Box<[F]>,
}

/// Framework-neutral graph matching `PyG` and `DGL` tensor conventions.
#[derive(Clone, Debug, PartialEq)]
pub struct Graph {
    /// Number of nodes.
    pub node_count: usize,
    /// Row-major `(2, edges)` signed indices: all sources, then all targets.
    pub edge_index: Box<[i64]>,
    /// Number of directed or undirected edge records.
    pub edge_count: usize,
    /// Dense `(nodes, features)` matrix.
    pub node_features: FeatureMatrix<NodeFeature>,
    /// Dense `(edges, features)` matrix.
    pub edge_features: FeatureMatrix<EdgeFeature>,
}

impl Graph {
    /// Graph export materialises framework-oriented indices and feature arrays.
    #[must_use]
    pub const fn cost(&self) -> crate::ExportCost {
        crate::ExportCost::Copy
    }
}

/// Invalid graph request or unavailable source data.
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum GraphError {
    /// A parameter is zero, non-finite or otherwise invalid.
    #[error("invalid graph construction parameter")]
    InvalidParameter,
    /// The requested feature is not meaningful for the selected node/edge mode.
    #[error("feature {feature} is unsupported for this graph configuration")]
    UnsupportedFeature {
        /// Stable feature label.
        feature: &'static str,
    },
    /// A required value is absent at one output row.
    #[error("feature {feature} is absent at row {row}")]
    MissingFeature {
        /// Stable feature label.
        feature: &'static str,
        /// Zero-based node or edge row.
        row: usize,
    },
    /// Periodic export was requested from a structure without a usable cell.
    #[error("periodic graph export requires a valid unit cell")]
    MissingCell,
    /// Topology rows and the selected coordinate block disagree in length.
    #[error("graph export requires one dense coordinate model")]
    RaggedCoordinates,
    /// Spatial indexing rejected coordinates, indices, cutoff or cell data.
    #[error(transparent)]
    Spatial(#[from] SpatialError),
    /// A graph index cannot be represented by the signed tensor convention.
    #[error("graph index exceeds signed 64-bit tensor limits")]
    IndexOverflow,
    /// Graph dimensions overflow or the host cannot reserve the requested buffer.
    #[error("graph tensor dimensions or allocation exceed host limits")]
    ResourceLimit,
}
