use pdbiox_core::contract::{AnalysisPolicy, PolicyField};

/// One completed analysis at one policy point.
#[derive(Clone, Debug, PartialEq)]
pub struct AuditRun<R> {
    /// The exact policy used.
    pub policy: AnalysisPolicy,
    /// The caller's unmodified analysis result.
    pub result: R,
}

/// One item that is not present under every audited policy.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SensitiveItem<I> {
    /// Stable identity supplied by the audited analysis.
    pub item: I,
    /// Zero-based run indices in which the item is present.
    pub present_in: Vec<usize>,
}

/// Sensitivity attributable to one varied policy dimension.
#[derive(Clone, Debug, PartialEq)]
pub struct DimensionSensitivity<I> {
    /// The field whose alternatives were compared.
    pub field: PolicyField,
    /// Items whose presence changes when this field changes while others stay fixed.
    pub sensitive_items: Vec<I>,
    /// Mean Jaccard loss across pairs differing only in this field.
    pub mean_change: f64,
}

/// Full result of a bounded policy audit.
#[derive(Clone, Debug, PartialEq)]
pub struct AuditReport<R, I> {
    /// Completed runs in the plan's deterministic order.
    pub runs: Vec<AuditRun<R>>,
    /// Fraction of distinct items present in every run; one for no items.
    pub stability: f64,
    /// Every non-invariant item and the runs that contain it.
    pub sensitive_items: Vec<SensitiveItem<I>>,
    /// Per-field sensitivity, in plan declaration order.
    pub dimensions: Vec<DimensionSensitivity<I>>,
}
