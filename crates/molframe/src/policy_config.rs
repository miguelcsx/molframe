//! Strict JSON and TOML loading for complete analysis policies.

use crate::document::{DocumentError, read_document};
use molframe_core::contract::{AnalysisPolicy, PolicyParseError};
use serde::Deserialize;
use std::path::{Path, PathBuf};

/// Failure while loading a declarative analysis policy.
#[derive(Debug, thiserror::Error)]
pub enum PolicyConfigError {
    /// Policy file could not be read.
    #[error("policy file I/O failed: {0}")]
    Io(#[from] std::io::Error),
    /// JSON syntax or schema was invalid.
    #[error("invalid JSON policy: {0}")]
    Json(#[from] serde_json::Error),
    /// TOML syntax or schema was invalid.
    #[error("invalid TOML policy: {0}")]
    Toml(#[from] toml::de::Error),
    /// File suffix does not name a supported configuration syntax.
    #[error("policy files must use a .json or .toml suffix")]
    UnsupportedFormat,
    /// A field used a value outside its documented vocabulary.
    #[error("invalid policy value for {field}: {value}")]
    InvalidValue {
        /// Canonical policy field name.
        field: &'static str,
        /// Rejected value.
        value: String,
    },
}

impl From<DocumentError> for PolicyConfigError {
    fn from(error: DocumentError) -> Self {
        match error {
            DocumentError::Io(error) => Self::Io(error),
            DocumentError::Json(error) => Self::Json(error),
            DocumentError::Toml(error) => Self::Toml(error),
            DocumentError::UnsupportedFormat => Self::UnsupportedFormat,
        }
    }
}

/// Strict top-level configuration shared by library-backed applications.
#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct ApplicationConfiguration {
    /// Analysis-policy overrides.
    pub policy: PolicyOverrides,
    /// Output preferences for user-facing applications.
    pub output: OutputConfiguration,
    /// Explicit chemistry resource locations.
    pub chem: ChemistryConfiguration,
}

/// Output preferences that an application may apply before explicit flags.
#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct OutputConfiguration {
    /// Stable output-format vocabulary interpreted by the application.
    pub format: Option<String>,
}

/// Explicit chemistry resources; merely loading configuration performs no I/O.
#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct ChemistryConfiguration {
    /// Location of a caller-managed CCD cache.
    pub ccd_cache: Option<PathBuf>,
    /// Exact release identifier for the configured CCD cache.
    pub ccd_version: Option<String>,
}

/// String vocabulary accepted by policy files and command-line overrides.
#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct PolicyOverrides {
    /// Assembly vocabulary: `asymmetric-unit`, `biological:ID`, `crystal:RADIUS`.
    pub assembly: Option<String>,
    /// Model vocabulary: `first`, `index:N`, `all`, `ensemble`.
    pub model: Option<String>,
    /// Alternate-conformation vocabulary.
    pub altloc: Option<String>,
    /// Identifier namespace vocabulary.
    pub identifiers: Option<String>,
    /// Missing-atom policy vocabulary.
    pub missing_atoms: Option<String>,
    /// Hydrogen policy vocabulary.
    pub hydrogens: Option<String>,
    /// Atom-equivalence vocabulary.
    pub atom_equivalence: Option<String>,
    /// Symmetry vocabulary.
    pub symmetry: Option<String>,
    /// Alignment vocabulary.
    pub alignment: Option<String>,
    /// Accumulation precision vocabulary.
    pub precision: Option<String>,
    /// Periodic-boundary vocabulary.
    pub periodic: Option<String>,
    /// van der Waals radii vocabulary.
    pub vdw_radii: Option<String>,
    /// Contact definition with its explicit numeric parameter.
    pub contact_def: Option<String>,
    /// Relative floating-point tolerance.
    pub float_tolerance_relative: Option<f64>,
    /// Absolute floating-point tolerance.
    pub float_tolerance_absolute: Option<f64>,
}

/// Loads and validates a complete policy document from JSON or TOML.
///
/// Omitted fields retain the stable named default. Unknown fields are errors.
///
/// # Errors
///
/// Returns [`PolicyConfigError`] for I/O, syntax, schema or value errors.
pub fn read_policy(path: impl AsRef<Path>) -> Result<AnalysisPolicy, PolicyConfigError> {
    read_configuration(path)?
        .policy
        .apply_to(AnalysisPolicy::default())
}

/// Loads a strict application configuration without accessing referenced resources.
///
/// # Errors
///
/// Returns [`PolicyConfigError`] for I/O, syntax or schema errors.
pub fn read_configuration(
    path: impl AsRef<Path>,
) -> Result<ApplicationConfiguration, PolicyConfigError> {
    Ok(read_document(path.as_ref())?)
}

impl PolicyOverrides {
    /// Applies only present values over an existing policy.
    ///
    /// # Errors
    ///
    /// Returns [`PolicyConfigError::InvalidValue`] for invalid vocabulary or
    /// numeric parameters.
    pub fn apply_to(self, mut policy: AnalysisPolicy) -> Result<AnalysisPolicy, PolicyConfigError> {
        apply_optional(&mut policy.assembly, self.assembly)?;
        apply_optional(&mut policy.model, self.model)?;
        apply_optional(&mut policy.altloc, self.altloc)?;
        apply_optional(&mut policy.identifiers, self.identifiers)?;
        apply_optional(&mut policy.missing_atoms, self.missing_atoms)?;
        apply_optional(&mut policy.hydrogens, self.hydrogens)?;
        apply_optional(&mut policy.atom_equivalence, self.atom_equivalence)?;
        apply_optional(&mut policy.symmetry, self.symmetry)?;
        apply_optional(&mut policy.alignment, self.alignment)?;
        apply_optional(&mut policy.precision, self.precision)?;
        apply_optional(&mut policy.periodic, self.periodic)?;
        apply_optional(&mut policy.vdw_radii, self.vdw_radii)?;
        apply_optional(&mut policy.contact_def, self.contact_def)?;
        if let Some(relative) = self.float_tolerance_relative {
            require_nonnegative_finite("float_tolerance_relative", relative)?;
            policy.float_tolerance.relative = relative;
        }
        if let Some(absolute) = self.float_tolerance_absolute {
            require_nonnegative_finite("float_tolerance_absolute", absolute)?;
            policy.float_tolerance.absolute = absolute;
        }
        Ok(policy)
    }
}

fn apply_optional<T>(target: &mut T, value: Option<String>) -> Result<(), PolicyConfigError>
where
    T: std::str::FromStr<Err = PolicyParseError>,
{
    if let Some(value) = value {
        *target = value.parse()?;
    }
    Ok(())
}

impl From<PolicyParseError> for PolicyConfigError {
    fn from(error: PolicyParseError) -> Self {
        Self::InvalidValue {
            field: error.field,
            value: error.value,
        }
    }
}

fn require_nonnegative_finite(field: &'static str, value: f64) -> Result<(), PolicyConfigError> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(PolicyConfigError::InvalidValue {
            field,
            value: value.to_string(),
        })
    }
}

#[cfg(test)]
#[path = "policy_config_tests.rs"]
mod tests;
