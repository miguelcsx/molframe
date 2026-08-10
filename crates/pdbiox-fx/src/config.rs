//! Strict TOML and JSON loading for declarative motif evaluation.

use crate::{
    AtomSite, Comparison, ComponentRole, ComponentSpec, Constraint, MissingVerdict, Motif,
    MotifError, NamedConstraint, VerdictProfile, VerdictRule,
};
use serde::Deserialize;
use std::collections::BTreeMap;
use std::path::Path;

/// A validated motif and the explicit profile that decides its measurements.
#[derive(Debug)]
pub struct EvaluationSpecification {
    /// Declarative chemistry and geometry.
    pub motif: Motif,
    /// Named, versioned decision rules.
    pub profile: VerdictProfile,
}

/// Failure while loading a functional-evaluation specification.
#[derive(Debug, thiserror::Error)]
pub enum SpecificationError {
    /// File access failed.
    #[error("functional specification I/O failed: {0}")]
    Io(#[from] std::io::Error),
    /// JSON syntax or schema failed.
    #[error("invalid JSON functional specification: {0}")]
    Json(#[from] serde_json::Error),
    /// TOML syntax or schema failed.
    #[error("invalid TOML functional specification: {0}")]
    Toml(#[from] toml::de::Error),
    /// Unsupported filename suffix.
    #[error("functional specifications must use a .toml or .json suffix")]
    UnsupportedFormat,
    /// An atom site was not written as component.atom.
    #[error("invalid atom site: {0}")]
    InvalidAtomSite(String),
    /// A component role was not recognised.
    #[error("invalid component role: {0}")]
    InvalidRole(String),
    /// A verdict comparison was malformed.
    #[error("invalid comparison for metric {0}")]
    InvalidComparison(String),
    /// The validated motif rejected the document.
    #[error("invalid motif: {0}")]
    Motif(#[from] MotifError),
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Document {
    components: BTreeMap<String, ComponentDocument>,
    constraints: Vec<ConstraintDocument>,
    profile: ProfileDocument,
}

#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
struct ComponentDocument {
    role: String,
    component_ids: Vec<String>,
    required_atoms: Vec<String>,
    equivalent_atoms: Vec<Vec<String>>,
}

#[derive(Debug, Deserialize)]
#[serde(tag = "kind", rename_all = "kebab-case", deny_unknown_fields)]
enum ConstraintDocument {
    Distance {
        name: String,
        first: String,
        second: String,
        target: f64,
        tolerance: f64,
    },
    Angle {
        name: String,
        atoms: [String; 3],
        target: f64,
        tolerance: f64,
    },
    Dihedral {
        name: String,
        atoms: [String; 4],
        target: f64,
        tolerance: f64,
    },
    Chirality {
        name: String,
        atoms: [String; 4],
        positive: bool,
    },
    Planarity {
        name: String,
        atoms: Vec<String>,
        tolerance: f64,
    },
    Coordination {
        name: String,
        centre: String,
        partners: Vec<String>,
        count: usize,
        max_distance: f64,
    },
    StericExclusion {
        name: String,
        first: String,
        second: String,
        min_distance: f64,
    },
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct ProfileDocument {
    id: String,
    missing: String,
    rules: Vec<RuleDocument>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RuleDocument {
    metric: String,
    comparison: String,
    values: Vec<f64>,
}

/// Loads and validates a functional-evaluation specification.
///
/// # Errors
///
/// Returns [`SpecificationError`] for I/O, syntax, schema or semantic errors.
pub fn read_evaluation_specification(
    path: impl AsRef<Path>,
) -> Result<EvaluationSpecification, SpecificationError> {
    let path = path.as_ref();
    let text = std::fs::read_to_string(path)?;
    let document: Document = match path.extension().and_then(std::ffi::OsStr::to_str) {
        Some(extension) if extension.eq_ignore_ascii_case("toml") => toml::from_str(&text)?,
        Some(extension) if extension.eq_ignore_ascii_case("json") => serde_json::from_str(&text)?,
        _ => return Err(SpecificationError::UnsupportedFormat),
    };
    document.validate()
}

impl Document {
    fn validate(self) -> Result<EvaluationSpecification, SpecificationError> {
        let components = self
            .components
            .into_iter()
            .map(|(name, document)| Ok((name.into_boxed_str(), document.validate()?)))
            .collect::<Result<Vec<_>, SpecificationError>>()?;
        let constraints = self
            .constraints
            .into_iter()
            .map(ConstraintDocument::validate)
            .collect::<Result<Vec<_>, _>>()?;
        let motif = Motif::new(components, constraints)?;
        let profile = self.profile.validate()?;
        Ok(EvaluationSpecification { motif, profile })
    }
}

impl ComponentDocument {
    fn validate(self) -> Result<ComponentSpec, SpecificationError> {
        let role = match self.role.as_str() {
            "residue" => ComponentRole::Residue,
            "ligand" => ComponentRole::Ligand,
            "cofactor" => ComponentRole::Cofactor,
            "metal" => ComponentRole::Metal,
            _ => return Err(SpecificationError::InvalidRole(self.role)),
        };
        let mut component = ComponentSpec::new(role);
        for id in self.component_ids {
            component = component.component(id);
        }
        for atom in self.required_atoms {
            component = component.require(atom);
        }
        for group in self.equivalent_atoms {
            component = component.equivalent(group);
        }
        Ok(component)
    }
}

impl ConstraintDocument {
    fn validate(self) -> Result<NamedConstraint, SpecificationError> {
        let constraint = match self {
            Self::Distance {
                name,
                first,
                second,
                target,
                tolerance,
            } => named(
                name,
                Constraint::Distance {
                    first: site(&first)?,
                    second: site(&second)?,
                    target,
                    tolerance,
                },
            ),
            Self::Angle {
                name,
                atoms,
                target,
                tolerance,
            } => named(
                name,
                Constraint::Angle {
                    atoms: sites(atoms)?,
                    target,
                    tolerance,
                },
            ),
            Self::Dihedral {
                name,
                atoms,
                target,
                tolerance,
            } => named(
                name,
                Constraint::Dihedral {
                    atoms: sites(atoms)?,
                    target,
                    tolerance,
                },
            ),
            Self::Chirality {
                name,
                atoms,
                positive,
            } => named(
                name,
                Constraint::Chirality {
                    atoms: sites(atoms)?,
                    positive,
                },
            ),
            Self::Planarity {
                name,
                atoms,
                tolerance,
            } => named(
                name,
                Constraint::Planarity {
                    atoms: site_list(&atoms)?,
                    tolerance,
                },
            ),
            Self::Coordination {
                name,
                centre,
                partners,
                count,
                max_distance,
            } => named(
                name,
                Constraint::Coordination {
                    centre: site(&centre)?,
                    partners: site_list(&partners)?,
                    count,
                    max_distance,
                },
            ),
            Self::StericExclusion {
                name,
                first,
                second,
                min_distance,
            } => named(
                name,
                Constraint::StericExclusion {
                    first: site(&first)?,
                    second: site(&second)?,
                    min_distance,
                },
            ),
        };
        Ok(constraint)
    }
}

fn named(name: String, constraint: Constraint) -> NamedConstraint {
    NamedConstraint {
        name: name.into(),
        constraint,
    }
}

fn site_list(values: &[String]) -> Result<Vec<AtomSite>, SpecificationError> {
    values.iter().map(|value| site(value)).collect()
}

impl ProfileDocument {
    fn validate(self) -> Result<VerdictProfile, SpecificationError> {
        if self.id.is_empty() || self.rules.is_empty() {
            return Err(SpecificationError::InvalidComparison(self.id));
        }
        let missing = match self.missing.as_str() {
            "indeterminate" => MissingVerdict::Indeterminate,
            "fail" => MissingVerdict::Fail,
            _ => return Err(SpecificationError::InvalidComparison("missing".to_owned())),
        };
        let rules = self
            .rules
            .into_iter()
            .map(RuleDocument::validate)
            .collect::<Result<Vec<_>, _>>()?;
        Ok(VerdictProfile::new(self.id, rules, missing))
    }
}

impl RuleDocument {
    fn validate(self) -> Result<VerdictRule, SpecificationError> {
        if self.metric.is_empty() || self.values.iter().any(|value| !value.is_finite()) {
            return Err(SpecificationError::InvalidComparison(self.metric));
        }
        let comparison = match (self.comparison.as_str(), self.values.as_slice()) {
            ("less-than", [value]) => Comparison::LessThan(*value),
            ("at-most", [value]) => Comparison::AtMost(*value),
            ("greater-than", [value]) => Comparison::GreaterThan(*value),
            ("at-least", [value]) => Comparison::AtLeast(*value),
            ("between", [low, high]) if low <= high => Comparison::Between(*low, *high),
            _ => return Err(SpecificationError::InvalidComparison(self.metric)),
        };
        Ok(VerdictRule {
            metric: self.metric.into(),
            comparison,
        })
    }
}

fn site(value: &str) -> Result<AtomSite, SpecificationError> {
    let Some((component, atom)) = value.split_once('.') else {
        return Err(SpecificationError::InvalidAtomSite(value.to_owned()));
    };
    if component.is_empty() || atom.is_empty() {
        return Err(SpecificationError::InvalidAtomSite(value.to_owned()));
    }
    Ok(AtomSite::new(component, atom))
}

fn sites<const N: usize>(values: [String; N]) -> Result<[AtomSite; N], SpecificationError> {
    let parsed = values.map(|value| site(&value));
    parsed
        .into_iter()
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .map_err(|_| SpecificationError::InvalidAtomSite("invalid site count".to_owned()))
}

#[cfg(test)]
#[path = "config_tests.rs"]
mod tests;
