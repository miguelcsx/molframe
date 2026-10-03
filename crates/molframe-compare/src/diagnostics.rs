//! Stable diagnostic codes for every error the comparison kernels return.

use crate::{
    CadConstructionError, CadError, CeError, CompareError, GovernedCompareError, MappedCompareError,
};
use molframe_core::{Code, Diagnostic, diagnostic_from};

diagnostic_from!(CompareError, |error| match error {
    CompareError::LengthMismatch { .. } | CompareError::WeightLengthMismatch { .. } => {
        Code::E5102
    }
    CompareError::Superpose(_) => Code::E5104,
    CompareError::MappingLimit { .. } => Code::E1901,
    CompareError::EmptyComponent | CompareError::NoComparablePairs => Code::E5003,
    CompareError::InvalidMapping => Code::E5102,
    CompareError::InvalidDistanceCutoff
    | CompareError::InvalidLengthScale
    | CompareError::InvalidScoreInput => Code::E5101,
    CompareError::Spatial(_) => Code::E4002,
    CompareError::UnsupportedNamespace => Code::E6103,
});

diagnostic_from!(CadError, |_error| Code::E5101);

diagnostic_from!(CadConstructionError, |error| match error {
    CadConstructionError::LengthMismatch => Code::E5102,
    CadConstructionError::Surface(_) => Code::E5101,
});

diagnostic_from!(CeError, |error| match error {
    CeError::InvalidOptions => Code::E5101,
    CeError::TooFewPoints { .. } => Code::E5103,
    CeError::MemoryLimit { .. } => Code::E7001,
    CeError::NoAlignment => Code::E5105,
    CeError::Superpose(_) => Code::E5104,
});

impl From<&GovernedCompareError> for Diagnostic {
    /// A governed failure names the code of the failure it wraps.
    fn from(error: &GovernedCompareError) -> Self {
        match error {
            GovernedCompareError::Compare(inner) => Self::from(inner),
            GovernedCompareError::Cad(inner) => Self::from(inner),
            GovernedCompareError::CadConstruction(inner) => Self::from(inner),
            GovernedCompareError::Ce(inner) => Self::from(inner),
            GovernedCompareError::Diagnostic(inner) => inner.clone(),
            GovernedCompareError::CoverageOverflow => {
                Self::new(Code::E1903).with_message(error.to_string())
            }
            GovernedCompareError::MissingData(_) => {
                Self::new(Code::E5102).with_message(error.to_string())
            }
            GovernedCompareError::UnsupportedPolicy(_) => {
                Self::new(Code::E6103).with_message(error.to_string())
            }
        }
    }
}

impl From<GovernedCompareError> for Diagnostic {
    fn from(error: GovernedCompareError) -> Self {
        Self::from(&error)
    }
}

impl From<&MappedCompareError> for Diagnostic {
    /// A mapping failure carries its own diagnostic; a score failure its kernel's.
    fn from(error: &MappedCompareError) -> Self {
        match error {
            MappedCompareError::Mapping(inner) => inner.clone(),
            MappedCompareError::Compare(inner) => Self::from(inner),
        }
    }
}

impl From<MappedCompareError> for Diagnostic {
    fn from(error: MappedCompareError) -> Self {
        Self::from(&error)
    }
}

#[cfg(test)]
#[path = "diagnostics_tests.rs"]
mod tests;
