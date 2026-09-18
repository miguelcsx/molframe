//! Strict JSON and TOML loading for complete analysis policies.

use molframe_core::contract::{
    AlignmentPolicy, AltlocPolicy, AnalysisPolicy, AssemblyChoice, ContactDefinition,
    EquivalencePolicy, HydrogenPolicy, MissingPolicy, ModelChoice, Namespace, PeriodicPolicy,
    Precision, RadiiSet, SymmetryPolicy,
};
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
    read_policy_overrides(path)?.apply_to(AnalysisPolicy::default())
}

/// Loads a partial policy document without applying a baseline.
///
/// # Errors
///
/// Returns [`PolicyConfigError`] for I/O, syntax or schema errors.
pub fn read_policy_overrides(path: impl AsRef<Path>) -> Result<PolicyOverrides, PolicyConfigError> {
    Ok(read_configuration(path)?.policy)
}

/// Loads a strict application configuration without accessing referenced resources.
///
/// # Errors
///
/// Returns [`PolicyConfigError`] for I/O, syntax or schema errors.
pub fn read_configuration(
    path: impl AsRef<Path>,
) -> Result<ApplicationConfiguration, PolicyConfigError> {
    let path = path.as_ref();
    let text = std::fs::read_to_string(path)?;
    let document: ApplicationConfiguration = match path
        .extension()
        .and_then(std::ffi::OsStr::to_str)
    {
        Some(extension) if extension.eq_ignore_ascii_case("json") => serde_json::from_str(&text)?,
        Some(extension) if extension.eq_ignore_ascii_case("toml") => toml::from_str(&text)?,
        _ => return Err(PolicyConfigError::UnsupportedFormat),
    };
    Ok(document)
}

impl PolicyOverrides {
    /// Applies only present values over an existing policy.
    ///
    /// # Errors
    ///
    /// Returns [`PolicyConfigError::InvalidValue`] for invalid vocabulary or
    /// numeric parameters.
    pub fn apply_to(self, mut policy: AnalysisPolicy) -> Result<AnalysisPolicy, PolicyConfigError> {
        apply_optional(&mut policy.assembly, self.assembly, parse_assembly)?;
        apply_optional(&mut policy.model, self.model, parse_model)?;
        apply_optional(&mut policy.altloc, self.altloc, parse_altloc)?;
        apply_optional(&mut policy.identifiers, self.identifiers, parse_identifiers)?;
        apply_optional(&mut policy.missing_atoms, self.missing_atoms, parse_missing)?;
        apply_optional(&mut policy.hydrogens, self.hydrogens, parse_hydrogens)?;
        apply_optional(
            &mut policy.atom_equivalence,
            self.atom_equivalence,
            parse_equivalence,
        )?;
        apply_optional(&mut policy.symmetry, self.symmetry, parse_symmetry)?;
        apply_optional(&mut policy.alignment, self.alignment, parse_alignment)?;
        apply_optional(&mut policy.precision, self.precision, parse_precision)?;
        apply_optional(&mut policy.periodic, self.periodic, parse_periodic)?;
        apply_optional(&mut policy.vdw_radii, self.vdw_radii, parse_radii)?;
        apply_optional(&mut policy.contact_def, self.contact_def, parse_contact)?;
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

fn apply_optional<T>(
    target: &mut T,
    value: Option<String>,
    parse: fn(&str) -> Result<T, PolicyConfigError>,
) -> Result<(), PolicyConfigError> {
    if let Some(value) = value {
        *target = parse(&value)?;
    }
    Ok(())
}

fn parse_assembly(value: &str) -> Result<AssemblyChoice, PolicyConfigError> {
    match value {
        "asymmetric-unit" => Ok(AssemblyChoice::AsymmetricUnit),
        _ if value.starts_with("biological:") => Ok(AssemblyChoice::Biological(
            required_suffix("assembly", value, "biological:")?.into(),
        )),
        _ if value.starts_with("crystal:") => Ok(AssemblyChoice::Crystal {
            radius: positive_f32("assembly", required_suffix("assembly", value, "crystal:")?)?,
        }),
        _ => invalid("assembly", value),
    }
}

fn parse_model(value: &str) -> Result<ModelChoice, PolicyConfigError> {
    match value {
        "first" => Ok(ModelChoice::First),
        "all" => Ok(ModelChoice::All),
        "ensemble" => Ok(ModelChoice::Ensemble),
        _ if value.starts_with("index:") => Ok(ModelChoice::Index(parse_u32(
            "model",
            required_suffix("model", value, "index:")?,
        )?)),
        _ => invalid("model", value),
    }
}

fn parse_altloc(value: &str) -> Result<AltlocPolicy, PolicyConfigError> {
    match value {
        "keep-all" => Ok(AltlocPolicy::KeepAll),
        "conformer-consistent" => Ok(AltlocPolicy::ConformerConsistent),
        "first" => Ok(AltlocPolicy::First),
        "highest-occupancy-per-residue" => Ok(AltlocPolicy::HighestOccupancyPerResidue),
        "highest-occupancy-per-atom" => Ok(AltlocPolicy::HighestOccupancyPerAtom),
        _ if value.starts_with("label:") => Ok(AltlocPolicy::Label(
            required_suffix("altloc", value, "label:")?.into(),
        )),
        _ => invalid("altloc", value),
    }
}

fn parse_identifiers(value: &str) -> Result<Namespace, PolicyConfigError> {
    parse_closed(
        "identifiers",
        value,
        &[
            ("label", Namespace::Label),
            ("auth", Namespace::Auth),
            ("explicit", Namespace::Explicit),
        ],
    )
}

fn parse_missing(value: &str) -> Result<MissingPolicy, PolicyConfigError> {
    parse_closed(
        "missing_atoms",
        value,
        &[
            ("ignore", MissingPolicy::Ignore),
            ("report", MissingPolicy::Report),
            ("indeterminate", MissingPolicy::Indeterminate),
            ("fail", MissingPolicy::Fail),
        ],
    )
}

fn parse_hydrogens(value: &str) -> Result<HydrogenPolicy, PolicyConfigError> {
    parse_closed(
        "hydrogens",
        value,
        &[
            ("explicit-only", HydrogenPolicy::ExplicitOnly),
            ("include-inferred", HydrogenPolicy::IncludeInferred),
            ("exclude", HydrogenPolicy::Exclude),
        ],
    )
}

fn parse_equivalence(value: &str) -> Result<EquivalencePolicy, PolicyConfigError> {
    parse_closed(
        "atom_equivalence",
        value,
        &[
            ("none", EquivalencePolicy::None),
            ("ccd", EquivalencePolicy::Ccd),
            ("explicit", EquivalencePolicy::Explicit),
        ],
    )
}

fn parse_symmetry(value: &str) -> Result<SymmetryPolicy, PolicyConfigError> {
    parse_closed(
        "symmetry",
        value,
        &[
            ("none", SymmetryPolicy::None),
            ("crystallographic", SymmetryPolicy::Crystallographic),
            ("biological-assembly", SymmetryPolicy::BiologicalAssembly),
        ],
    )
}

fn parse_alignment(value: &str) -> Result<AlignmentPolicy, PolicyConfigError> {
    match value {
        "none" => Ok(AlignmentPolicy::None),
        "global" => Ok(AlignmentPolicy::Global),
        "local" => Ok(AlignmentPolicy::Local),
        _ if value.starts_with("explicit:") => Ok(AlignmentPolicy::Explicit(
            required_suffix("alignment", value, "explicit:")?.into(),
        )),
        _ => invalid("alignment", value),
    }
}

fn parse_precision(value: &str) -> Result<Precision, PolicyConfigError> {
    parse_closed(
        "precision",
        value,
        &[("f32", Precision::F32), ("f64", Precision::F64)],
    )
}

fn parse_periodic(value: &str) -> Result<PeriodicPolicy, PolicyConfigError> {
    parse_closed(
        "periodic",
        value,
        &[
            ("none", PeriodicPolicy::None),
            ("pbc", PeriodicPolicy::Pbc),
            ("minimum-image", PeriodicPolicy::MinimumImage),
        ],
    )
}

fn parse_radii(value: &str) -> Result<RadiiSet, PolicyConfigError> {
    parse_closed(
        "vdw_radii",
        value,
        &[
            ("bondi", RadiiSet::Bondi),
            ("amber-united", RadiiSet::AmberUnited),
            ("charmm", RadiiSet::Charmm),
            ("alvarez", RadiiSet::Alvarez),
        ],
    )
}

fn parse_contact(value: &str) -> Result<ContactDefinition, PolicyConfigError> {
    if let Some(tolerance) = value.strip_prefix("distance:") {
        return Ok(ContactDefinition::DistanceCutoff {
            tolerance: nonnegative_f32("contact_def", tolerance)?,
        });
    }
    if let Some(probe) = value.strip_prefix("surface:") {
        return Ok(ContactDefinition::SurfaceBased {
            probe: positive_f32("contact_def", probe)?,
        });
    }
    invalid("contact_def", value)
}

fn parse_closed<T: Clone>(
    field: &'static str,
    value: &str,
    choices: &[(&str, T)],
) -> Result<T, PolicyConfigError> {
    choices
        .iter()
        .find(|(name, _)| *name == value)
        .map(|(_, result)| result.clone())
        .ok_or_else(|| invalid_error(field, value))
}

fn required_suffix<'a>(
    field: &'static str,
    value: &'a str,
    prefix: &str,
) -> Result<&'a str, PolicyConfigError> {
    value
        .strip_prefix(prefix)
        .filter(|suffix| !suffix.is_empty())
        .ok_or_else(|| invalid_error(field, value))
}

fn parse_u32(field: &'static str, value: &str) -> Result<u32, PolicyConfigError> {
    value.parse().map_err(|_| invalid_error(field, value))
}

fn positive_f32(field: &'static str, value: &str) -> Result<f32, PolicyConfigError> {
    let parsed = nonnegative_f32(field, value)?;
    if parsed > 0.0 {
        Ok(parsed)
    } else {
        invalid(field, value)
    }
}

fn nonnegative_f32(field: &'static str, value: &str) -> Result<f32, PolicyConfigError> {
    let parsed: f32 = value.parse().map_err(|_| invalid_error(field, value))?;
    if parsed.is_finite() && parsed >= 0.0 {
        Ok(parsed)
    } else {
        invalid(field, value)
    }
}

fn require_nonnegative_finite(field: &'static str, value: f64) -> Result<(), PolicyConfigError> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        invalid(field, &value.to_string())
    }
}

fn invalid<T>(field: &'static str, value: &str) -> Result<T, PolicyConfigError> {
    Err(invalid_error(field, value))
}

fn invalid_error(field: &'static str, value: &str) -> PolicyConfigError {
    PolicyConfigError::InvalidValue {
        field,
        value: value.to_owned(),
    }
}

#[cfg(test)]
#[path = "policy_config_tests.rs"]
mod tests;
