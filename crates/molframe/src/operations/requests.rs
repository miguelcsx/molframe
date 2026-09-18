//! Requests that carry their own validation.
//!
//! Most families build their request from a plain descriptor and validate it in
//! the executor. These three cannot: a selection query has to be compiled and a
//! contact cutoff checked before a plan will accept them, so each constructor
//! here is the only way to obtain one and every instance is valid by
//! construction.

use molframe_core::contract::AnalysisPolicy;
use molframe_query::Query;
use molframe_spatial::SpatialBackend;

use super::plan::value::ExecutionPlanError;

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

/// One coordinate-array slot in a plan's input, typed so a request cannot
/// silently confuse a slot number with an atom or model index.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CoordinateSlot(usize);

impl CoordinateSlot {
    /// Creates a slot handle for the caller's plan input.
    #[must_use]
    pub const fn new(slot: usize) -> Self {
        Self(slot)
    }

    /// The slot number in the plan input's array list.
    #[must_use]
    pub const fn slot(self) -> usize {
        self.0
    }
}

/// A typed request for a coordinate-array RMSD calculation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RmsdRequest {
    mobile: CoordinateSlot,
    reference: CoordinateSlot,
}

impl RmsdRequest {
    /// Creates a request over two coordinate-array slots.
    #[must_use]
    pub const fn new(mobile: CoordinateSlot, reference: CoordinateSlot) -> Self {
        Self { mobile, reference }
    }

    /// Mobile coordinate-array slot.
    #[must_use]
    pub const fn mobile(&self) -> CoordinateSlot {
        self.mobile
    }

    /// Reference coordinate-array slot.
    #[must_use]
    pub const fn reference(&self) -> CoordinateSlot {
        self.reference
    }
}

impl From<RmsdRequest> for super::plan::value::PlanOperation {
    fn from(value: RmsdRequest) -> Self {
        Self::Rmsd(value)
    }
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
        let query = Query::compile(source.as_ref())
            .map_err(|findings| ExecutionPlanError::Selection(findings.into()))?;
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

impl From<SelectionRequest> for super::plan::value::PlanOperation {
    fn from(value: SelectionRequest) -> Self {
        Self::Selection(value)
    }
}

impl From<ContactsRequest> for super::plan::value::PlanOperation {
    fn from(value: ContactsRequest) -> Self {
        Self::Contacts(Box::new(value))
    }
}
