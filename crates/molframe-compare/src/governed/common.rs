//! Shared governed comparison errors and coverage construction.

use crate::{CadConstructionError, CadError, CeError, CompareError};
use molframe_core::Diagnostic;
use molframe_core::contract::{
    AlgorithmId, Analysis, AnalysisPolicy, Coverage, MissingPolicy, ParameterValue, Status,
};

/// A governed comparison failure.
#[derive(Debug, thiserror::Error)]
pub enum GovernedCompareError {
    /// The scientific comparison kernel refused its inputs.
    #[error(transparent)]
    Compare(#[from] CompareError),
    /// Contact-area scoring rejected an invalid area map.
    #[error(transparent)]
    Cad(#[from] CadError),
    /// Surface contact-area construction refused its inputs.
    #[error(transparent)]
    CadConstruction(#[from] CadConstructionError),
    /// Combinatorial Extension rejected its inputs or found no alignment.
    #[error(transparent)]
    Ce(#[from] CeError),
    /// Chemistry-backed sequence or identifier resolution failed.
    #[error("comparison input resolution failed: {0}")]
    Diagnostic(Diagnostic),
    /// The public coverage counters cannot represent the input cardinality.
    #[error("comparison input count exceeds u32 coverage limits")]
    CoverageOverflow,
    /// Policy rejects an incomplete explicit correspondence.
    #[error("comparison correspondence omits {0} intended inputs")]
    MissingData(u32),
    /// The requested policy cannot be represented by this explicit workflow.
    #[error("unsupported comparison policy field: {0}")]
    UnsupportedPolicy(&'static str),
}

pub(crate) fn float(value: impl Into<f64>) -> ParameterValue {
    ParameterValue::Float(value.into().to_bits())
}

pub(crate) fn complete<T>(
    value: T,
    count: usize,
    policy: &AnalysisPolicy,
    name: &'static str,
) -> Result<Analysis<T>, GovernedCompareError> {
    let intended = u32::try_from(count).map_err(|_| GovernedCompareError::CoverageOverflow)?;
    let mut analysis = Analysis::complete(value, Coverage::complete(intended), policy);
    analysis.provenance = analysis
        .provenance
        .with_algorithm(AlgorithmId::new(name, "1"));
    Ok(analysis)
}

pub(crate) fn covered<T>(
    value: T,
    intended: usize,
    used: usize,
    policy: &AnalysisPolicy,
    name: &'static str,
) -> Result<Analysis<T>, GovernedCompareError> {
    let intended = u32::try_from(intended).map_err(|_| GovernedCompareError::CoverageOverflow)?;
    let used = u32::try_from(used).map_err(|_| GovernedCompareError::CoverageOverflow)?;
    let missing = intended
        .checked_sub(used)
        .ok_or(GovernedCompareError::CoverageOverflow)?;
    let status = match (missing, policy.missing_atoms) {
        (0, _) => Status::Complete,
        (_, MissingPolicy::Ignore | MissingPolicy::Report) => Status::Partial,
        (_, MissingPolicy::Indeterminate) => Status::Indeterminate,
        (_, MissingPolicy::Fail) => return Err(GovernedCompareError::MissingData(missing)),
        _ => return Err(GovernedCompareError::UnsupportedPolicy("missing_atoms")),
    };
    let mut analysis = Analysis::complete(
        value,
        Coverage {
            intended,
            used,
            missing,
            ambiguous: 0,
        },
        policy,
    );
    analysis.status = status;
    analysis.provenance = analysis
        .provenance
        .with_algorithm(AlgorithmId::new(name, "1"));
    Ok(analysis)
}
