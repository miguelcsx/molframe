//! The decision behind an ambiguous status.

use molframe_core::contract::{
    AltlocPolicy, AnalysisPolicy, Assumption, AssumptionSource, Impact, PolicyField,
};

/// Names the altloc policy that produced an ambiguous result, so the status is
/// never left without a stated cause.
pub(super) fn ambiguity_assumption(
    policy: &AnalysisPolicy,
    analysed: u32,
    input: u32,
    left_out: u32,
) -> Assumption {
    let source = if matches!(policy.altloc, AltlocPolicy::ConformerConsistent) {
        AssumptionSource::ProfileDefault
    } else {
        AssumptionSource::Explicit
    };
    Assumption::new(
        PolicyField::Altloc,
        format!(
            "alternate conformations are present and the altloc policy `{:?}` is applied: {analysed} of {input} atoms analysed, {left_out} left out as other conformations",
            policy.altloc,
        ),
        source,
        Impact::Unmeasured,
    )
}
