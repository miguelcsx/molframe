//! Stable diagnostic codes for the chemistry parsers and queries.

use crate::{AutomorphismLimit, Mol2Error, MolError, SmartsDataError, SmartsError};
use molframe_core::{Code, diagnostic_from};

diagnostic_from!(MolError, |error| match error {
    MolError::Malformed => Code::E1201,
    MolError::MetadataLength => Code::E5102,
    MolError::BondEndpoint => Code::E3006,
    MolError::Unrepresentable(_) => Code::E4105,
});

diagnostic_from!(Mol2Error, |error| match error {
    Mol2Error::Malformed => Code::E1201,
    Mol2Error::CountMismatch | Mol2Error::MetadataLength => Code::E5102,
    Mol2Error::InvalidIdentifier => Code::E2002,
    Mol2Error::Unrepresentable => Code::E4105,
});

diagnostic_from!(SmartsError, |_error| Code::E1301, |diagnostic, error| {
    diagnostic.with_context("position", error.position.to_string())
});

diagnostic_from!(SmartsDataError, |_error| Code::E4003);

diagnostic_from!(AutomorphismLimit, |_error| Code::E1901);

#[cfg(test)]
#[path = "diagnostics_tests.rs"]
mod tests;
