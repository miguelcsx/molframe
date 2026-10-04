//! Stable diagnostic codes for every error the validation kernels return.

use crate::policy_execution::{
    AltlocOccupancyKernelError, BFactorKernelError, CcdCompletenessKernelError, GovernedMapError,
    LigandGeometryKernelError, PlaneRestraintKernelError,
};
use crate::{
    AltlocOccupancyError, BFactorError, CompletenessError, NucleicGeometryError, PlanarityError,
    RamachandranError, RealSpaceCorrelationError, ReferenceError, RotamerError,
};
use molframe_core::{Code, Diagnostic, diagnostic_from};

diagnostic_from!(AltlocOccupancyError, |error| match error {
    AltlocOccupancyError::InvalidOptions => Code::E5101,
    AltlocOccupancyError::ExplicitNamespace
    | AltlocOccupancyError::UnsupportedNamespace
    | AltlocOccupancyError::MissingAtomName => Code::E6103,
});
diagnostic_from!(CompletenessError, |_error| Code::E6103);
diagnostic_from!(PlanarityError, |error| match error {
    PlanarityError::InvalidTolerance => Code::E5101,
    PlanarityError::Geometry(inner) => Diagnostic::from(inner).code(),
});
diagnostic_from!(ReferenceError, |_error| Code::E5101);
diagnostic_from!(RamachandranError, |error| match error {
    RamachandranError::CoverageOverflow => Code::E1903,
    RamachandranError::NoBasins
    | RamachandranError::DuplicateBasin
    | RamachandranError::DuplicateDistribution
    | RamachandranError::OutlierBasin
    | RamachandranError::EmptyDistribution
    | RamachandranError::Probability => Code::E5101,
    RamachandranError::Reference(inner) => Diagnostic::from(inner).code(),
    RamachandranError::Roles(inner) => inner.code(),
});
diagnostic_from!(RotamerError, |error| match error {
    RotamerError::CoverageOverflow => Code::E1903,
    RotamerError::InvalidProfile
    | RotamerError::InvalidOptions
    | RotamerError::InvalidPath { .. } => Code::E5101,
    RotamerError::Provider(inner) => inner.code(),
    RotamerError::Reference(inner) => Diagnostic::from(inner).code(),
});
diagnostic_from!(RealSpaceCorrelationError, |error| match error {
    RealSpaceCorrelationError::IncompatibleGrid => Code::E5102,
    RealSpaceCorrelationError::InsufficientSamples => Code::E5103,
    RealSpaceCorrelationError::NonFiniteDensity => Code::E5101,
    RealSpaceCorrelationError::ZeroVariance => Code::E5104,
});
diagnostic_from!(NucleicGeometryError, |error| match error {
    NucleicGeometryError::InvalidPolicy => Code::E5101,
    NucleicGeometryError::Torsions(inner) => Diagnostic::from(inner).code(),
    NucleicGeometryError::Geometry(inner) => Diagnostic::from(inner).code(),
    NucleicGeometryError::Provider(inner) => inner.code(),
    NucleicGeometryError::MissingComponentIdentity { .. }
    | NucleicGeometryError::MissingComponent { .. } => Code::E6103,
    NucleicGeometryError::AmbiguousRole { .. } => Code::E5105,
});
diagnostic_from!(BFactorError, |error| match error {
    BFactorError::InvalidThreshold | BFactorError::InvalidTls => Code::E5101,
    BFactorError::NoObservations => Code::E5103,
});
diagnostic_from!(PlaneRestraintKernelError, |error| match error {
    PlaneRestraintKernelError::Validation(inner) => Diagnostic::from(inner).code(),
    PlaneRestraintKernelError::CoverageOverflow => Code::E1903,
});
diagnostic_from!(LigandGeometryKernelError, |error| match error {
    LigandGeometryKernelError::CoverageOverflow => Code::E1903,
});
diagnostic_from!(AltlocOccupancyKernelError, |error| match error {
    AltlocOccupancyKernelError::Validation(inner) => Diagnostic::from(inner).code(),
    AltlocOccupancyKernelError::CoverageOverflow => Code::E1903,
});
diagnostic_from!(CcdCompletenessKernelError, |error| match error {
    CcdCompletenessKernelError::Diagnostic(inner) => inner.code(),
    CcdCompletenessKernelError::CoverageOverflow => Code::E1903,
});
diagnostic_from!(BFactorKernelError, |error| match error {
    BFactorKernelError::Validation(inner) => Diagnostic::from(inner).code(),
    BFactorKernelError::CoverageOverflow => Code::E1903,
});
diagnostic_from!(GovernedMapError, |error| match error {
    GovernedMapError::Correlation(inner) => Diagnostic::from(inner).code(),
    GovernedMapError::CoverageOverflow => Code::E1903,
    GovernedMapError::MissingSamples(_) => Code::E5103,
    GovernedMapError::UnsupportedPolicy(_) => Code::E6103,
});

#[cfg(test)]
#[path = "diagnostics_tests.rs"]
mod tests;
