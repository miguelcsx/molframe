//! The words a policy decision is spelled with.
//!
//! One vocabulary serves every surface that names a decision in text: a
//! configuration file, a command-line flag and the Python keyword arguments all
//! parse and print through these implementations, so a spelling cannot mean two
//! things. Names are kebab-case; an underscore is accepted for a hyphen, because
//! Python spells them that way, and a decision that carries a value is written
//! `name:value` (`biological:1`, `crystal:12.5`, `index:3`, `label:A`).

use super::{
    AlignmentPolicy, AltlocPolicy, AssemblyChoice, ContactDefinition, EquivalencePolicy,
    HydrogenPolicy, MissingPolicy, ModelChoice, Namespace, PeriodicPolicy, Precision, RadiiSet,
    SymmetryPolicy,
};
use std::fmt;
use std::str::FromStr;

/// A named choice (a policy decision, a spatial backend, a file format) was
/// spelled with a word outside its vocabulary.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[error("invalid policy value for {field}: {value} (expected {expected})")]
pub struct PolicyParseError {
    /// Which decision was being read.
    pub field: &'static str,
    /// What was written.
    pub value: String,
    /// What would have been accepted, for the message.
    pub expected: String,
}

crate::diagnostic_from!(
    PolicyParseError,
    |_error| crate::Code::E6104,
    |diagnostic, error| diagnostic
        .about_field(error.field)
        .with_context("value", error.value.clone())
        .with_context("expected", error.expected.clone())
);

impl PolicyParseError {
    /// A refusal of `value` for the choice named `field`, which would have
    /// accepted `expected`.
    #[must_use]
    pub fn new(field: &'static str, value: &str, expected: &str) -> Self {
        Self {
            field,
            value: value.to_owned(),
            expected: expected.to_owned(),
        }
    }
}

fn refused(field: &'static str, value: &str, expected: &str) -> PolicyParseError {
    PolicyParseError::new(field, value, expected)
}

/// `value` with underscores read as hyphens, the spelling every vocabulary
/// compares against.
#[must_use]
pub fn canonical_spelling(value: &str) -> String {
    value.replace('_', "-")
}

fn canonical(value: &str) -> String {
    canonical_spelling(value)
}

/// The payload of a `prefix:payload` spelling, which must not be empty.
fn payload<'a>(
    field: &'static str,
    value: &'a str,
    prefix: &str,
) -> Result<Option<&'a str>, PolicyParseError> {
    match value.split_once(':') {
        Some((head, rest)) if canonical(head) == prefix => {
            if rest.is_empty() {
                Err(refused(field, value, &format!("a value after '{prefix}:'")))
            } else {
                Ok(Some(rest))
            }
        }
        _ => Ok(None),
    }
}

fn finite_float(
    field: &'static str,
    value: &str,
    text: &str,
    allow_zero: bool,
) -> Result<f32, PolicyParseError> {
    let expected = if allow_zero {
        "a finite non-negative number"
    } else {
        "a finite positive number"
    };
    let parsed: f32 = text.parse().map_err(|_| refused(field, value, expected))?;
    if parsed.is_finite() && (parsed > 0.0 || (allow_zero && parsed == 0.0)) {
        Ok(parsed)
    } else {
        Err(refused(field, value, expected))
    }
}

/// Declares the spelling of a decision that is a plain closed set of words.
macro_rules! closed_vocabulary {
    ($type:ty, $field:literal, [$(($variant:path, $name:literal)),+ $(,)?]) => {
        impl $type {
            /// The words this decision can be spelled with, in declaration order.
            pub const NAMES: &'static [&'static str] = &[$($name),+];

            /// The canonical spelling of this choice.
            #[must_use]
            pub const fn name(&self) -> &'static str {
                match self {
                    $($variant => $name,)+
                    #[allow(unreachable_patterns)]
                    _ => "custom",
                }
            }
        }

        impl fmt::Display for $type {
            fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                formatter.write_str(self.name())
            }
        }

        impl FromStr for $type {
            type Err = PolicyParseError;

            fn from_str(value: &str) -> Result<Self, Self::Err> {
                match canonical(value).as_str() {
                    $($name => Ok($variant),)+
                    _ => Err(refused($field, value, &Self::NAMES.join(", "))),
                }
            }
        }
    };
}

closed_vocabulary!(
    Namespace,
    "identifiers",
    [
        (Namespace::Label, "label"),
        (Namespace::Auth, "auth"),
        (Namespace::Explicit, "explicit"),
    ]
);
closed_vocabulary!(
    MissingPolicy,
    "missing_atoms",
    [
        (MissingPolicy::Ignore, "ignore"),
        (MissingPolicy::Report, "report"),
        (MissingPolicy::Indeterminate, "indeterminate"),
        (MissingPolicy::Fail, "fail"),
    ]
);
closed_vocabulary!(
    HydrogenPolicy,
    "hydrogens",
    [
        (HydrogenPolicy::ExplicitOnly, "explicit-only"),
        (HydrogenPolicy::IncludeInferred, "include-inferred"),
        (HydrogenPolicy::Exclude, "exclude"),
    ]
);
closed_vocabulary!(
    EquivalencePolicy,
    "atom_equivalence",
    [
        (EquivalencePolicy::None, "none"),
        (EquivalencePolicy::Ccd, "ccd"),
        (EquivalencePolicy::Explicit, "explicit"),
    ]
);
closed_vocabulary!(
    SymmetryPolicy,
    "symmetry",
    [
        (SymmetryPolicy::None, "none"),
        (SymmetryPolicy::Crystallographic, "crystallographic"),
        (SymmetryPolicy::BiologicalAssembly, "biological-assembly"),
    ]
);
closed_vocabulary!(
    Precision,
    "precision",
    [(Precision::F32, "f32"), (Precision::F64, "f64"),]
);
closed_vocabulary!(
    PeriodicPolicy,
    "periodic",
    [
        (PeriodicPolicy::None, "none"),
        (PeriodicPolicy::Pbc, "pbc"),
        (PeriodicPolicy::MinimumImage, "minimum-image"),
    ]
);
closed_vocabulary!(
    RadiiSet,
    "vdw_radii",
    [
        (RadiiSet::Bondi, "bondi"),
        (RadiiSet::AmberUnited, "amber-united"),
        (RadiiSet::Charmm, "charmm"),
        (RadiiSet::Alvarez, "alvarez"),
    ]
);

impl FromStr for AltlocPolicy {
    type Err = PolicyParseError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        if let Some(label) = payload("altloc", value, "label")? {
            return Ok(Self::Label(label.into()));
        }
        match canonical(value).as_str() {
            "keep-all" => Ok(Self::KeepAll),
            "conformer-consistent" => Ok(Self::ConformerConsistent),
            "first" => Ok(Self::First),
            "highest-occupancy-per-residue" => Ok(Self::HighestOccupancyPerResidue),
            "highest-occupancy-per-atom" => Ok(Self::HighestOccupancyPerAtom),
            _ => Err(refused(
                "altloc",
                value,
                "keep-all, conformer-consistent, first, highest-occupancy-per-residue, \
                 highest-occupancy-per-atom or label:<id>",
            )),
        }
    }
}

impl fmt::Display for AltlocPolicy {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::KeepAll => formatter.write_str("keep-all"),
            Self::ConformerConsistent => formatter.write_str("conformer-consistent"),
            Self::First => formatter.write_str("first"),
            Self::HighestOccupancyPerResidue => {
                formatter.write_str("highest-occupancy-per-residue")
            }
            Self::HighestOccupancyPerAtom => formatter.write_str("highest-occupancy-per-atom"),
            Self::Label(label) => write!(formatter, "label:{label}"),
        }
    }
}

impl FromStr for AssemblyChoice {
    type Err = PolicyParseError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        if let Some(id) = payload("assembly", value, "biological")? {
            return Ok(Self::Biological(id.into()));
        }
        if let Some(radius) = payload("assembly", value, "crystal")? {
            return Ok(Self::Crystal {
                radius: finite_float("assembly", value, radius, false)?,
            });
        }
        match canonical(value).as_str() {
            "asymmetric-unit" => Ok(Self::AsymmetricUnit),
            _ => Err(refused(
                "assembly",
                value,
                "asymmetric-unit, biological:<id> or crystal:<radius>",
            )),
        }
    }
}

impl fmt::Display for AssemblyChoice {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::AsymmetricUnit => formatter.write_str("asymmetric-unit"),
            Self::Biological(id) => write!(formatter, "biological:{id}"),
            Self::Crystal { radius } => write!(formatter, "crystal:{radius}"),
        }
    }
}

impl FromStr for ModelChoice {
    type Err = PolicyParseError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        if let Some(index) = payload("model", value, "index")? {
            return index
                .parse()
                .map(Self::Index)
                .map_err(|_| refused("model", value, "index:<non-negative integer>"));
        }
        match canonical(value).as_str() {
            "first" => Ok(Self::First),
            "all" => Ok(Self::All),
            "ensemble" => Ok(Self::Ensemble),
            _ => Err(refused("model", value, "first, all, ensemble or index:<n>")),
        }
    }
}

impl fmt::Display for ModelChoice {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::First => formatter.write_str("first"),
            Self::Index(index) => write!(formatter, "index:{index}"),
            Self::All => formatter.write_str("all"),
            Self::Ensemble => formatter.write_str("ensemble"),
        }
    }
}

impl FromStr for AlignmentPolicy {
    type Err = PolicyParseError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        if let Some(selection) = payload("alignment", value, "explicit")? {
            return Ok(Self::Explicit(selection.into()));
        }
        match canonical(value).as_str() {
            "none" => Ok(Self::None),
            "global" => Ok(Self::Global),
            "local" => Ok(Self::Local),
            _ => Err(refused(
                "alignment",
                value,
                "none, global, local or explicit:<selection>",
            )),
        }
    }
}

impl fmt::Display for AlignmentPolicy {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::None => formatter.write_str("none"),
            Self::Explicit(selection) => write!(formatter, "explicit:{selection}"),
            Self::Global => formatter.write_str("global"),
            Self::Local => formatter.write_str("local"),
        }
    }
}

impl FromStr for ContactDefinition {
    type Err = PolicyParseError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        if let Some(tolerance) = payload("contact_def", value, "distance")? {
            return Ok(Self::DistanceCutoff {
                tolerance: finite_float("contact_def", value, tolerance, true)?,
            });
        }
        if let Some(probe) = payload("contact_def", value, "surface")? {
            return Ok(Self::SurfaceBased {
                probe: finite_float("contact_def", value, probe, false)?,
            });
        }
        Err(refused(
            "contact_def",
            value,
            "distance:<tolerance> or surface:<probe>",
        ))
    }
}

impl fmt::Display for ContactDefinition {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DistanceCutoff { tolerance } => write!(formatter, "distance:{tolerance}"),
            Self::SurfaceBased { probe } => write!(formatter, "surface:{probe}"),
        }
    }
}

#[cfg(test)]
#[path = "vocabulary_tests.rs"]
mod tests;
