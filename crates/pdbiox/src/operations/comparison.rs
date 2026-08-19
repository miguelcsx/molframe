//! Typed coordinate comparison requests used by the facade plan.

use super::PlanOperation;

/// A comparison kernel that consumes two coordinate arrays.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ComparisonMetric {
    /// Local distance-difference test with a caller-selected inclusion radius.
    Lddt {
        /// Largest reference distance included in the local domain.
        inclusion_radius: f64,
    },
    /// Length-normalised fitted score.
    TmScore,
    /// Global distance test over the standard TS cutoffs.
    GdtTs,
    /// Global distance test over the high-accuracy cutoffs.
    GdtHa,
}

/// A validated coordinate comparison request.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ComparisonRequest {
    mobile: usize,
    reference: usize,
    metric: ComparisonMetric,
}

impl ComparisonRequest {
    /// Creates a comparison request for two coordinate-array slots.
    ///
    /// # Panics
    ///
    /// This constructor does not validate the slot numbers because slots are
    /// owned by the caller's plan input. lDDT's radius is validated before
    /// execution and reports a native comparison error.
    #[must_use]
    pub const fn new(mobile: usize, reference: usize, metric: ComparisonMetric) -> Self {
        Self {
            mobile,
            reference,
            metric,
        }
    }

    /// Mobile coordinate-array slot.
    #[must_use]
    pub const fn mobile(&self) -> usize {
        self.mobile
    }

    /// Reference coordinate-array slot.
    #[must_use]
    pub const fn reference(&self) -> usize {
        self.reference
    }

    /// Selected comparison metric.
    #[must_use]
    pub const fn metric(&self) -> ComparisonMetric {
        self.metric
    }
}

/// The typed score returned by a comparison plan node.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ComparisonResult {
    /// Metric used to produce the score.
    pub metric: ComparisonMetric,
    /// Native score value.
    pub value: f64,
}

/// Executes one comparison without materialising either input array.
pub(crate) fn execute(
    request: &ComparisonRequest,
    mobile: &[[f32; 3]],
    reference: &[[f32; 3]],
) -> Result<ComparisonResult, pdbiox_compare::CompareError> {
    let value = match request.metric {
        ComparisonMetric::Lddt { inclusion_radius } => pdbiox_compare::lddt_with_options(
            mobile,
            reference,
            &pdbiox_compare::LddtOptions::standard(inclusion_radius),
        )?,
        ComparisonMetric::TmScore => pdbiox_compare::tm_score(mobile, reference)?,
        ComparisonMetric::GdtTs => pdbiox_compare::gdt_ts(mobile, reference)?,
        ComparisonMetric::GdtHa => pdbiox_compare::gdt_ha(mobile, reference)?,
    };
    Ok(ComparisonResult {
        metric: request.metric,
        value,
    })
}

impl From<ComparisonRequest> for PlanOperation {
    fn from(value: ComparisonRequest) -> Self {
        Self::Comparison(value)
    }
}
