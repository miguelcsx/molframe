//! Stable diagnostic codes for every error the analysis kernels return.
//!
//! Wrapped errors name the code of the failure they wrap. The rest name the
//! kind of thing that went wrong: an invalid parameter, inputs that do not
//! correspond, a method that did not converge, a missing prerequisite or a
//! resource limit; the message carries which.

use crate::policy_execution::{
    GovernedAnalysisError, GovernedNativeError, PhysicalKernelError, StandaloneAnalysisError,
};
use crate::{
    AnmError, BasePairError, CationPiError, DsspBinaryError, DsspError, FragmentMappingError,
    GnmError, HelicalError, HseError, HydrogenBondError, NativeError, NmdError,
    NucleicTorsionError, PiStackingError, PolymerError, PoreError, PotentialError, RadialError,
};
use molframe_core::{Code, Diagnostic, diagnostic_from};

diagnostic_from!(HelicalError, |error| match error {
    HelicalError::InvalidTolerance => Code::E5101,
    HelicalError::InvalidFrame => Code::E5104,
});
diagnostic_from!(NativeError, |error| match error {
    NativeError::AtomCountMismatch { .. } => Code::E5102,
    NativeError::Spatial(inner) => Diagnostic::from(inner).code(),
});
diagnostic_from!(PiStackingError, |error| match error {
    PiStackingError::InvalidOptions => Code::E5101,
    PiStackingError::Geometry(inner) => Diagnostic::from(inner).code(),
});
diagnostic_from!(CationPiError, |error| match error {
    CationPiError::InvalidOptions => Code::E5101,
    CationPiError::RingGeometry(inner) => Diagnostic::from(inner).code(),
});
diagnostic_from!(RadialError, |error| match error {
    RadialError::InvalidBounds
    | RadialError::NoBins
    | RadialError::InvalidVolume
    | RadialError::InvalidGroups => Code::E5101,
    RadialError::PairCountOverflow => Code::E1903,
    RadialError::Spatial(inner) => Diagnostic::from(inner).code(),
});
diagnostic_from!(AnmError, |error| match error {
    AnmError::InvalidOptions => Code::E5101,
    AnmError::TooFewSites | AnmError::InsufficientModes { .. } => Code::E5103,
    AnmError::MemoryLimit { .. } => Code::E7001,
    AnmError::Convergence { .. } | AnmError::CoincidentSites { .. } => Code::E5104,
    AnmError::Spatial(inner) => Diagnostic::from(inner).code(),
});
diagnostic_from!(GnmError, |error| match error {
    GnmError::InvalidOptions => Code::E5101,
    GnmError::TooFewSites | GnmError::InsufficientModes { .. } => Code::E5103,
    GnmError::MemoryLimit { .. } => Code::E7001,
    GnmError::Convergence { .. } => Code::E5104,
    GnmError::Spatial(inner) => Diagnostic::from(inner).code(),
});
diagnostic_from!(NmdError, |error| match error {
    NmdError::NotUtf8 | NmdError::InvalidNumber { .. } | NmdError::ModeBeforeCoordinates => {
        Code::E7101
    }
    NmdError::UnalignedVectors { .. } | NmdError::SiteMismatch { .. } => Code::E5102,
});
diagnostic_from!(HydrogenBondError, |error| match error {
    HydrogenBondError::InvalidOptions => Code::E5101,
    HydrogenBondError::MissingChemistry => Code::E6103,
    HydrogenBondError::MissingCell => Code::E5004,
    HydrogenBondError::Spatial(inner) => Diagnostic::from(inner).code(),
});
diagnostic_from!(DsspError, |error| match error {
    DsspError::MissingRoleAnnotation => Code::E6103,
    DsspError::InvalidOptions => Code::E5101,
    DsspError::AmbiguousRole { .. } => Code::E5105,
});
diagnostic_from!(PotentialError, |error| match error {
    PotentialError::MisalignedInput => Code::E5102,
    PotentialError::NonFiniteInput | PotentialError::InvalidGrid => Code::E5101,
    PotentialError::GridTooLarge => Code::E7001,
    PotentialError::Spatial(inner) => Diagnostic::from(inner).code(),
    PotentialError::WorkerPanicked => Code::E9001,
});
diagnostic_from!(HseError, |error| match error {
    HseError::MissingRoleAnnotation => Code::E6103,
    HseError::AmbiguousRole { .. } => Code::E5105,
    HseError::Spatial(inner) => Diagnostic::from(inner).code(),
});
diagnostic_from!(NucleicTorsionError, |error| match error {
    NucleicTorsionError::MissingRoles => Code::E6103,
    NucleicTorsionError::AmbiguousRole { .. } => Code::E5105,
});
diagnostic_from!(PoreError, |error| match error {
    PoreError::LengthMismatch { .. } => Code::E5102,
    PoreError::EmptyAtoms => Code::E5103,
    PoreError::InvalidAtoms | PoreError::InvalidOptions => Code::E5101,
    PoreError::GridTooLarge | PoreError::MemoryLimit { .. } => Code::E7001,
});
diagnostic_from!(PolymerError, |error| match error {
    PolymerError::TooShort => Code::E5103,
    PolymerError::NonFinite => Code::E5101,
    PolymerError::CoincidentSites => Code::E5104,
});
diagnostic_from!(DsspBinaryError, |error| match error {
    DsspBinaryError::Io(_) | DsspBinaryError::Exit { .. } => Code::E7101,
    DsspBinaryError::Parse(findings) => match findings.first() {
        Some(first) => first.code(),
        None => Code::E7101,
    },
});
diagnostic_from!(BasePairError, |error| match error {
    BasePairError::InvalidOptions => Code::E5101,
    BasePairError::MissingComponent { .. } => Code::E6103,
    BasePairError::Provider(inner) => inner.code(),
    BasePairError::HydrogenBond(inner) => Diagnostic::from(inner).code(),
});
diagnostic_from!(FragmentMappingError, |error| match error {
    FragmentMappingError::InvalidLibrary | FragmentMappingError::InvalidInput => Code::E5101,
    FragmentMappingError::Superpose(inner) => Diagnostic::from(inner).code(),
});

/// The diagnostic for a governed-execution failure, with `kernel` classifying
/// the kernel's own error type.
///
/// A function rather than a blanket `From`, because a generic wrapper that
/// requires `From` of its payload makes every conversion in this crate recurse
/// through it.
pub fn governed_diagnostic<E: std::fmt::Display>(
    error: &GovernedAnalysisError<E>,
    kernel: impl FnOnce(&E) -> Diagnostic,
) -> Diagnostic {
    let code = match error {
        GovernedAnalysisError::Trajectory(inner) => Diagnostic::from(inner).code(),
        // The structure's own finding says what is wrong with it; the wrapper
        // would only say that something is.
        GovernedAnalysisError::InvalidStructure(findings) => match findings.first() {
            Some(first) => return first.clone(),
            None => Code::E9001,
        },
        GovernedAnalysisError::Kernel(inner) => return kernel(inner),
        GovernedAnalysisError::MissingData { .. } => Code::E5103,
        GovernedAnalysisError::CoverageOverflow => Code::E1903,
        GovernedAnalysisError::MultipleModelsRequested
        | GovernedAnalysisError::UnsupportedPolicyValue(_) => Code::E6103,
        GovernedAnalysisError::MissingFrameOutput
        | GovernedAnalysisError::InvalidCoverage { .. } => Code::E9001,
    };
    Diagnostic::new(code).with_message(error.to_string())
}

impl From<&GovernedNativeError> for Diagnostic {
    fn from(error: &GovernedNativeError) -> Self {
        let code = match error {
            GovernedNativeError::Memory(_) => Code::E7001,
            GovernedNativeError::Native(inner) => Self::from(inner).code(),
            GovernedNativeError::InvalidTarget(findings) => match findings.first() {
                Some(first) => first.code(),
                None => Code::E9001,
            },
            GovernedNativeError::SourceIndexOverflow(_) => Code::E1903,
            GovernedNativeError::TargetAtomOutOfRange { .. } => Code::E5102,
            GovernedNativeError::MissingTargetModel(_) => Code::E6009,
            GovernedNativeError::MultipleModelsRequested
            | GovernedNativeError::PeriodicUnsupported
            | GovernedNativeError::UnsupportedPolicy(_) => Code::E6103,
        };
        Self::new(code).with_message(error.to_string())
    }
}

/// As [`governed_diagnostic`], for a physical kernel's error.
pub fn physical_diagnostic<E: std::fmt::Display>(
    error: &PhysicalKernelError<E>,
    kernel: impl FnOnce(&E) -> Diagnostic,
) -> Diagnostic {
    let code = match error {
        PhysicalKernelError::Kernel(inner) => return kernel(inner),
        PhysicalKernelError::MissingCell => Code::E5004,
        PhysicalKernelError::AllImagesUnsupported => Code::E6103,
        PhysicalKernelError::SideInputLength => Code::E5102,
    };
    Diagnostic::new(code).with_message(error.to_string())
}

/// As [`governed_diagnostic`], for a standalone analysis's error.
pub fn standalone_diagnostic<E: std::fmt::Display>(
    error: &StandaloneAnalysisError<E>,
    kernel: impl FnOnce(&E) -> Diagnostic,
) -> Diagnostic {
    let code = match error {
        StandaloneAnalysisError::Kernel(inner) => return kernel(inner),
        StandaloneAnalysisError::CoverageOverflow => Code::E1903,
        StandaloneAnalysisError::MissingData { .. } => Code::E5103,
        StandaloneAnalysisError::UnsupportedMissingPolicy => Code::E6103,
    };
    Diagnostic::new(code).with_message(error.to_string())
}

#[cfg(test)]
#[path = "diagnostics_tests.rs"]
mod tests;
