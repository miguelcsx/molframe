//! The typed vocabulary a plan is written in: what it can hold, what it
//! produces, and how it fails.
//!
//! Each operation variant is one family's request, and each value variant is
//! that family's typed result. The two enums are the whole cross-family
//! contract — a family adds a variant here and its kernel, and the executor in
//! [`Plan`](crate::Plan) is the only place that has to know both.

use molframe_analysis::Contact;
use molframe_core::contract::Analysis;
use molframe_core::diagnostic::{Diagnostic, Findings};
use molframe_geom::SuperposeError;
use std::fmt;
use thiserror::Error;

use super::super::requests::{ContactsRequest, RmsdRequest, SelectionRequest};

/// One operation accepted by [`Plan`](crate::Plan).
#[derive(Clone, Debug)]
pub enum PlanOperation {
    /// Compiled query selection over a structure snapshot.
    Selection(SelectionRequest),
    /// Selection-scoped atom-contact operation.
    Contacts(Box<ContactsRequest>),
    /// RMSD over two borrowed coordinate-array slots.
    Rmsd(RmsdRequest),
    /// A policy-bound structure analysis or validation operation.
    Structure(Box<super::super::StructureRequest>),
    /// A policy-bound physical analysis operation.
    Physical(Box<super::super::PhysicalRequest>),
    /// A typed coordinate or ensemble geometry operation.
    Geometry(Box<super::super::geometry::GeometryRequest>),
    /// A fixed-radius spatial search over borrowed coordinates.
    Spatial(Box<super::super::spatial::SpatialRequest>),
    /// A surface-area calculation over borrowed coordinates and atom data.
    #[cfg(feature = "surface")]
    Surface(Box<super::super::surface::SurfaceRequest>),
    /// Explicit distance-based covalent-bond inference.
    BondInference(crate::BondInference),
    /// A comparison score over two borrowed coordinate-array slots.
    #[cfg(feature = "compare")]
    Comparison(super::super::ComparisonRequest),
    /// A governed trajectory analysis over one borrowed frame-array slot.
    #[cfg(feature = "traj")]
    Trajectory(Box<super::super::trajectory::TrajectoryRequest>),
}

// Every family's request converts into an operation, so a plan is written in
// one voice: `plan.add(id, request)`.
impl From<super::super::StructureRequest> for PlanOperation {
    fn from(value: super::super::StructureRequest) -> Self {
        Self::Structure(Box::new(value))
    }
}

impl From<super::super::PhysicalRequest> for PlanOperation {
    fn from(value: super::super::PhysicalRequest) -> Self {
        Self::Physical(Box::new(value))
    }
}

impl From<super::super::geometry::GeometryRequest> for PlanOperation {
    fn from(value: super::super::geometry::GeometryRequest) -> Self {
        Self::Geometry(Box::new(value))
    }
}

impl From<super::super::spatial::SpatialRequest> for PlanOperation {
    fn from(value: super::super::spatial::SpatialRequest) -> Self {
        Self::Spatial(Box::new(value))
    }
}

#[cfg(feature = "surface")]
impl From<super::super::surface::SurfaceRequest> for PlanOperation {
    fn from(value: super::super::surface::SurfaceRequest) -> Self {
        Self::Surface(Box::new(value))
    }
}

impl From<crate::BondInference> for PlanOperation {
    fn from(value: crate::BondInference) -> Self {
        Self::BondInference(value)
    }
}

/// A typed result value produced by one plan node.
#[derive(Clone, Debug)]
pub enum PlanValue {
    /// Query evaluation with the selected atom set and diagnostics.
    Selection(Box<crate::Evaluation>),
    /// Atom contacts with the full analysis contract.
    Contacts(Box<Analysis<Vec<Contact>>>),
    /// Coordinate RMSD.
    Rmsd(f64),
    /// A policy-bound structure result with typed scientific output.
    Structure(Box<super::super::StructureValue>),
    /// A typed physical-analysis result.
    Physical(Box<super::super::PhysicalValue>),
    /// A typed coordinate or ensemble geometry result.
    Geometry(Box<super::super::geometry::GeometryValue>),
    /// A typed spatial-search result.
    Spatial(Box<super::super::spatial::SpatialValue>),
    /// A typed solvent-accessible or buried-surface result.
    #[cfg(feature = "surface")]
    Surface(Box<super::super::surface::SurfaceValue>),
    /// A structure snapshot with inferred connectivity and skipped atoms.
    BondInference(Box<crate::BondInferenceReport>),
    /// A typed coordinate comparison score.
    #[cfg(feature = "compare")]
    Comparison(super::super::ComparisonResult),
    /// A typed governed trajectory analysis result.
    #[cfg(feature = "traj")]
    Trajectory(Box<super::super::trajectory::TrajectoryValue>),
}

/// One named plan result.
#[derive(Clone, Debug)]
pub struct PlanResultEntry {
    /// Stable user-provided operation identifier.
    pub id: Box<str>,
    /// Typed native result.
    pub value: PlanValue,
}

/// Results from one deterministic plan execution.
#[derive(Clone, Debug)]
pub struct PlanResult {
    /// Results in sorted operation-id order.
    pub entries: Vec<PlanResultEntry>,
    /// Number of reusable spatial indices retained during execution.
    pub cached_index_count: usize,
}

/// Errors raised while building or executing a facade plan.
#[derive(Debug, Error)]
pub enum ExecutionPlanError {
    /// An operation id was repeated.
    #[error("duplicate operation id: {0}")]
    DuplicateId(Box<str>),
    /// A required input was not supplied.
    #[error("operation {operation:?} requires {input}")]
    MissingInput {
        /// Operation id.
        operation: Box<str>,
        /// Missing input class.
        input: &'static str,
    },
    /// An array slot did not exist.
    #[error("operation {operation:?} references array slot {slot}")]
    ArraySlot {
        /// Operation id.
        operation: Box<str>,
        /// Referenced array slot.
        slot: usize,
    },
    /// A scalar slot did not exist.
    #[error("operation {operation:?} references scalar slot {slot}")]
    ScalarSlot {
        /// Operation id.
        operation: Box<str>,
        /// Referenced scalar slot.
        slot: usize,
    },
    /// A frame slot did not exist.
    #[error("operation {operation:?} references frame slot {slot}")]
    FrameSlot {
        /// Operation id.
        operation: Box<str>,
        /// Referenced frame slot.
        slot: usize,
    },
    /// An atom-index slot did not exist.
    #[error("operation {operation:?} references index slot {slot}")]
    IndexSlot {
        /// Operation id.
        operation: Box<str>,
        /// Referenced atom-index slot.
        slot: usize,
    },
    /// A request was invalid before execution.
    #[error("invalid operation request: {0}")]
    InvalidRequest(Box<str>),
    /// Spatial planning or execution failed.
    #[error("spatial execution failed: {0:?}")]
    Spatial(Diagnostic),
    /// A fixed-radius spatial kernel rejected its input.
    #[error("spatial kernel failed: {0}")]
    SpatialKernel(Box<str>),
    /// A borrowed floating-point input did not exist.
    #[error("operation {operation:?} references float slot {slot}")]
    FloatSlot {
        /// Operation id.
        operation: Box<str>,
        /// Referenced floating-point slot.
        slot: usize,
    },
    /// A borrowed boolean input did not exist.
    #[cfg(feature = "surface")]
    #[error("operation {operation:?} references mask slot {slot}")]
    MaskSlot {
        /// Operation id.
        operation: Box<str>,
        /// Referenced boolean slot.
        slot: usize,
    },
    /// A surface kernel rejected its input.
    #[cfg(feature = "surface")]
    #[error("surface kernel failed: {0}")]
    SurfaceKernel(Box<str>),
    /// Chemistry planning or execution failed.
    #[error("chemistry execution failed: {0}")]
    Chemistry(Findings),
    /// Query selection failed with one or more diagnostics.
    #[error("selection execution failed: {0}")]
    Selection(Findings),
    /// A governed kernel failed.
    #[error("governed operation failed: {0}")]
    Governed(Box<str>),
    /// RMSD failed in the geometry kernel.
    #[error("RMSD failed: {0:?}")]
    Rmsd(SuperposeError),
    /// A dense geometry matrix could not be allocated.
    #[error("geometry matrix failed: {0}")]
    Matrix(molframe_geom::MatrixError),
    /// A small-matrix eigensolver failed.
    #[error("geometry eigendecomposition failed: {0:?}")]
    Eigen(molframe_geom::EigenError),
    /// Ensemble fluctuation input was invalid.
    #[error("geometry fluctuation failed: {0:?}")]
    Fluctuation(molframe_geom::FluctuationError),
    /// A comparison score failed in its native kernel.
    #[cfg(feature = "compare")]
    #[error("comparison failed: {0}")]
    Comparison(#[from] molframe_compare::CompareError),
    /// Borrowed trajectory-frame dimensions were invalid.
    #[cfg(feature = "traj")]
    #[error("trajectory frame input failed: {0}")]
    TrajectoryInput(#[from] molframe_traj::FrameViewError),
    /// A governed trajectory analysis failed.
    #[cfg(feature = "traj")]
    #[error(transparent)]
    Trajectory(#[from] molframe_traj::GovernedEnsembleError),
}

impl fmt::Display for PlanOperation {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Selection(_) => formatter.write_str("selection"),
            Self::Contacts(_) => formatter.write_str("contacts"),
            Self::Rmsd(_) => formatter.write_str("rmsd"),
            Self::Structure(_) => formatter.write_str("structure"),
            Self::Physical(_) => formatter.write_str("physical"),
            Self::Geometry(_) => formatter.write_str("geometry"),
            Self::Spatial(_) => formatter.write_str("spatial"),
            #[cfg(feature = "surface")]
            Self::Surface(_) => formatter.write_str("surface"),
            Self::BondInference(_) => formatter.write_str("bond_inference"),
            #[cfg(feature = "compare")]
            Self::Comparison(_) => formatter.write_str("comparison"),
            #[cfg(feature = "traj")]
            Self::Trajectory(_) => formatter.write_str("trajectory"),
        }
    }
}
