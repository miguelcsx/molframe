//! Bounded, deterministic sensitivity audits over [`AnalysisPolicy`] values.
//!
//! A caller declares the defensible alternatives for each policy field, inspects
//! the resulting cost, and then supplies an analysis plus a projection of its
//! result into individually identifiable items. The audit reports both global
//! stability and which items depend on each varied decision.

#![forbid(unsafe_code)]

mod batch;
mod engine;
mod numeric;
mod plan;
mod report;
mod value;

pub use batch::{BatchAudit, BatchDimension, audit_batch};
pub use engine::audit;
pub use plan::{AuditPlan, PlanError, PolicyDimension, PolicySpace};
pub use report::{AuditReport, AuditRun, DimensionSensitivity, SensitiveItem};
pub use value::PolicyValue;
