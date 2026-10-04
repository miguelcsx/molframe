//! Bounded, deterministic sensitivity audits over
//! [`AnalysisPolicy`](molframe_core::contract::AnalysisPolicy) values.
//!
//! A caller declares the defensible alternatives for each policy field, inspects
//! the resulting cost, and supplies an analysis. The audit runs every combination
//! and asks how far the answer moved, by a distance suited to what the answer is:
//! overlap of sets, error of a number, difference of vectors or rankings or
//! networks, or whether a conclusion flipped ([`OutcomeMetric`]).
//!
//! The variation among runs is then split among the decisions and their
//! interactions ([`decompose`]), so a decision that matters only together with
//! another is not mistaken for one that does not matter. A plan that varies a
//! field the analysis never applies is refused, and a combination of decisions
//! that contradict each other cannot be planned.

#![forbid(unsafe_code)]

mod batch;
mod design;
mod engine;
mod metric;
mod numeric;
mod outcomes;
mod plan;
mod report;
mod value;

pub use batch::{BatchAudit, BatchDimension, audit_batch};
pub use design::{Decomposition, Interaction, MainEffect, decompose};
pub use engine::audit;
pub use metric::{
    CategoricalFlip, Graph, GraphDifference, GraphPart, OutcomeMetric, RankingDistance,
    ScalarError, ScalarMode, SetOverlap, VectorDifference, VectorMode,
};
pub use outcomes::{AnalysisAudit, AuditError, OutcomeAudit, audit_analyses, audit_outcomes};
pub use plan::{AuditPlan, PlanError, PolicyDimension, PolicySpace};
pub use report::{AuditReport, AuditRun, DimensionSensitivity, SensitiveItem};
pub use value::PolicyValue;
