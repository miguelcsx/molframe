//! Owned operation requests and typed plan result values.

use pdbiox_analysis::Contact;
use pdbiox_core::contract::{Analysis, AnalysisPolicy};
use pdbiox_core::diagnostic::Diagnostic;
use pdbiox_core::structure::Structure;
use pdbiox_geom::SuperposeError;
use pdbiox_query::Query;
use pdbiox_spatial::SpatialBackend;
use std::fmt;
use thiserror::Error;

/// A typed request for a selection-scoped atom-contact calculation.
#[derive(Clone, Debug)]
pub struct ContactsRequest {
    left: Box<str>,
    right: Box<str>,
    left_query: Query,
    right_query: Query,
    cutoff: f32,
    backend: SpatialBackend,
    policy: AnalysisPolicy,
}

impl ContactsRequest {
    /// Creates a validated atom-contact request from two selection queries.
    ///
    /// # Errors
    ///
    /// Returns [`ExecutionPlanError::InvalidRequest`] for empty or invalid
    /// selections, or an invalid cutoff.
    pub fn new(
        left: impl Into<Box<str>>,
        right: impl Into<Box<str>>,
        cutoff: f32,
        backend: SpatialBackend,
        policy: AnalysisPolicy,
    ) -> Result<Self, ExecutionPlanError> {
        let left = left.into();
        let right = right.into();
        if left.is_empty() || right.is_empty() {
            return Err(ExecutionPlanError::InvalidRequest(
                "contact selections must not be empty".into(),
            ));
        }
        if !cutoff.is_finite() || cutoff <= 0.0 {
            return Err(ExecutionPlanError::InvalidRequest(
                "cutoff must be finite and positive".into(),
            ));
        }
        let left_query = Query::compile(&left).map_err(|findings| {
            ExecutionPlanError::InvalidRequest(format!("left selection: {findings:?}").into())
        })?;
        let right_query = Query::compile(&right).map_err(|findings| {
            ExecutionPlanError::InvalidRequest(format!("right selection: {findings:?}").into())
        })?;
        Ok(Self {
            left,
            right,
            left_query,
            right_query,
            cutoff,
            backend,
            policy,
        })
    }

    /// Left selection query.
    #[must_use]
    pub fn left(&self) -> &str {
        &self.left
    }

    /// Right selection query.
    #[must_use]
    pub fn right(&self) -> &str {
        &self.right
    }

    /// Compiled left selection query.
    #[must_use]
    pub fn left_query(&self) -> &Query {
        &self.left_query
    }

    /// Compiled right selection query.
    #[must_use]
    pub fn right_query(&self) -> &Query {
        &self.right_query
    }

    /// Distance cutoff in Angstroms.
    #[must_use]
    pub const fn cutoff(&self) -> f32 {
        self.cutoff
    }

    /// Requested spatial backend.
    #[must_use]
    pub const fn backend(&self) -> SpatialBackend {
        self.backend
    }

    /// Analysis policy retained by the request.
    #[must_use]
    pub fn policy(&self) -> &AnalysisPolicy {
        &self.policy
    }
}

/// A typed request for a coordinate-array RMSD calculation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RmsdRequest {
    mobile: usize,
    reference: usize,
}

/// A compiled structure-selection request evaluated by the native query engine.
#[derive(Clone, Debug)]
pub struct SelectionRequest {
    query: Query,
    policy: AnalysisPolicy,
}

impl SelectionRequest {
    /// Compiles textual query syntax and retains the policy for execution.
    ///
    /// # Errors
    ///
    /// Returns [`ExecutionPlanError::Selection`] when the query source cannot
    /// be compiled.
    pub fn new(
        source: impl AsRef<str>,
        policy: AnalysisPolicy,
    ) -> Result<Self, ExecutionPlanError> {
        let query = Query::compile(source.as_ref()).map_err(ExecutionPlanError::Selection)?;
        Ok(Self { query, policy })
    }

    /// Wraps an already compiled query without reparsing it.
    #[must_use]
    pub fn from_query(query: Query, policy: AnalysisPolicy) -> Self {
        Self { query, policy }
    }

    /// Compiled query retained by this request.
    #[must_use]
    pub fn query(&self) -> &Query {
        &self.query
    }

    /// Analysis policy retained by this request.
    #[must_use]
    pub fn policy(&self) -> &AnalysisPolicy {
        &self.policy
    }
}

impl RmsdRequest {
    /// Creates an array-slot request.
    #[must_use]
    pub const fn new(mobile: usize, reference: usize) -> Self {
        Self { mobile, reference }
    }

    /// Mobile array slot.
    #[must_use]
    pub const fn mobile(&self) -> usize {
        self.mobile
    }

    /// Reference array slot.
    #[must_use]
    pub const fn reference(&self) -> usize {
        self.reference
    }
}

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
    Structure(Box<super::StructureRequest>),
    /// A policy-bound physical analysis operation.
    Physical(Box<super::PhysicalRequest>),
    /// A typed coordinate or ensemble geometry operation.
    Geometry(Box<super::geometry::GeometryRequest>),
    /// A fixed-radius spatial search over borrowed coordinates.
    Spatial(Box<super::spatial::SpatialRequest>),
    /// A surface-area calculation over borrowed coordinates and atom data.
    #[cfg(feature = "surface")]
    Surface(Box<super::surface::SurfaceRequest>),
    /// Explicit distance-based covalent-bond inference.
    BondInference(crate::BondInference),
    /// A comparison score over two borrowed coordinate-array slots.
    #[cfg(feature = "compare")]
    Comparison(super::ComparisonRequest),
    /// A governed trajectory analysis over one borrowed frame-array slot.
    #[cfg(feature = "traj")]
    Trajectory(Box<super::trajectory::TrajectoryRequest>),
}

/// One borrowed coordinate array supplied to a plan execution.
#[derive(Clone, Copy, Debug)]
pub struct CoordinateInput<'a> {
    /// C-contiguous coordinate rows. The facade never materialises this slice.
    pub positions: &'a [[f32; 3]],
}

/// A borrowed flattened view over a C-contiguous `(frames, atoms, 3)` array.
#[derive(Clone, Copy, Debug)]
pub struct FrameInput<'a> {
    /// Flattened coordinate rows, retained by the binding owner.
    pub positions: &'a [[f32; 3]],
    /// Number of trajectory frames represented by `positions`.
    pub frame_count: usize,
    /// Number of atoms in each frame.
    pub atom_count: usize,
}

/// A borrowed C-contiguous atom-index array supplied to a plan execution.
#[derive(Clone, Copy, Debug)]
pub struct IndexInput<'a> {
    /// Zero-based atom indices. The facade never materialises this slice.
    pub indices: &'a [usize],
}

/// A borrowed one-dimensional `f64` input such as per-atom masses.
#[derive(Clone, Copy, Debug)]
pub struct ScalarInput<'a> {
    /// C-contiguous scalar values.
    pub values: &'a [f64],
}

/// A borrowed one-dimensional `f32` input used by atom-aligned analyses.
#[derive(Clone, Copy, Debug)]
pub struct FloatInput<'a> {
    /// C-contiguous floating-point values.
    pub values: &'a [f32],
}

/// A borrowed one-dimensional boolean input used by surface operations.
#[cfg(feature = "surface")]
#[derive(Clone, Copy, Debug)]
pub struct MaskInput<'a> {
    /// C-contiguous boolean values.
    pub values: &'a [bool],
}

/// Inputs shared by a plan execution.
#[derive(Clone, Copy, Debug, Default)]
pub struct PlanInput<'a> {
    /// Structure snapshot used by structure-bound operations.
    pub structure: Option<&'a Structure>,
    /// Borrowed arrays used by array-bound operations.
    pub arrays: &'a [CoordinateInput<'a>],
    /// Borrowed scalar arrays used by weighted geometry operations.
    pub scalars: &'a [ScalarInput<'a>],
    /// Borrowed floating-point arrays used by atom-aligned analyses.
    pub floats: &'a [FloatInput<'a>],
    /// Borrowed boolean arrays used by surface operations.
    #[cfg(feature = "surface")]
    pub masks: &'a [MaskInput<'a>],
    /// Borrowed frame arrays used by ensemble geometry operations.
    pub frames: &'a [FrameInput<'a>],
    /// Borrowed atom-index arrays used by trajectory analyses.
    pub indices: &'a [IndexInput<'a>],
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
    Structure(Box<super::StructureValue>),
    /// A typed physical-analysis result.
    Physical(Box<super::PhysicalValue>),
    /// A typed coordinate or ensemble geometry result.
    Geometry(Box<super::geometry::GeometryValue>),
    /// A typed spatial-search result.
    Spatial(Box<super::spatial::SpatialValue>),
    /// A typed solvent-accessible or buried-surface result.
    #[cfg(feature = "surface")]
    Surface(Box<super::surface::SurfaceValue>),
    /// A structure snapshot with inferred connectivity and skipped atoms.
    BondInference(Box<crate::BondInferenceReport>),
    /// A typed coordinate comparison score.
    #[cfg(feature = "compare")]
    Comparison(super::ComparisonResult),
    /// A typed governed trajectory analysis result.
    #[cfg(feature = "traj")]
    Trajectory(Box<super::trajectory::TrajectoryValue>),
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
    #[error("chemistry execution failed: {0:?}")]
    Chemistry(Diagnostic),
    /// Query selection failed with one or more diagnostics.
    #[error("selection execution failed: {0:?}")]
    Selection(Vec<Diagnostic>),
    /// A governed kernel failed.
    #[error("governed operation failed: {0}")]
    Governed(Box<str>),
    /// RMSD failed in the geometry kernel.
    #[error("RMSD failed: {0:?}")]
    Rmsd(SuperposeError),
    /// A dense geometry matrix could not be allocated.
    #[error("geometry matrix failed: {0}")]
    Matrix(pdbiox_geom::MatrixError),
    /// A small-matrix eigensolver failed.
    #[error("geometry eigendecomposition failed: {0:?}")]
    Eigen(pdbiox_geom::EigenError),
    /// Ensemble fluctuation input was invalid.
    #[error("geometry fluctuation failed: {0:?}")]
    Fluctuation(pdbiox_geom::FluctuationError),
    /// A comparison score failed in its native kernel.
    #[cfg(feature = "compare")]
    #[error("comparison failed: {0}")]
    Comparison(#[from] pdbiox_compare::CompareError),
    /// Borrowed trajectory-frame dimensions were invalid.
    #[cfg(feature = "traj")]
    #[error("trajectory frame input failed: {0}")]
    TrajectoryInput(#[from] pdbiox_traj::FrameViewError),
    /// A governed trajectory analysis failed.
    #[cfg(feature = "traj")]
    #[error(transparent)]
    Trajectory(#[from] pdbiox_traj::GovernedEnsembleError),
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
