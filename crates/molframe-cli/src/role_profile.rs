//! Strict caller-authored polymer role configuration; no atom names are implicit.

use molframe::chemistry::{ComponentKind, PolymerAtomRole, PolymerRoleProfile, PolymerRoleRule};
use serde::{Deserialize, Deserializer};
use std::path::Path;

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Document {
    profile_id: Box<str>,
    rules: Vec<Rule>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Rule {
    atom_name: Box<str>,
    #[serde(deserialize_with = "role")]
    role: PolymerAtomRole,
    component_id: Option<Box<str>>,
    #[serde(default, deserialize_with = "component_kind")]
    component_kind: Option<ComponentKind>,
}

fn role<'de, D: Deserializer<'de>>(deserializer: D) -> Result<PolymerAtomRole, D::Error> {
    let code = i64::deserialize(deserializer)?;
    PolymerAtomRole::from_code(code)
        .ok_or_else(|| serde::de::Error::custom("invalid polymer role bits"))
}

fn component_kind<'de, D: Deserializer<'de>>(
    deserializer: D,
) -> Result<Option<ComponentKind>, D::Error> {
    Option::<i64>::deserialize(deserializer)?
        .map(|code| {
            ComponentKind::from_code(code)
                .ok_or_else(|| serde::de::Error::custom("invalid component kind code"))
        })
        .transpose()
}

pub(crate) fn read(path: &Path) -> Result<PolymerRoleProfile, Box<dyn std::error::Error>> {
    let text = std::fs::read_to_string(path)?;
    let document: Document = match path.extension().and_then(std::ffi::OsStr::to_str) {
        Some(extension) if extension.eq_ignore_ascii_case("json") => serde_json::from_str(&text)?,
        Some(extension) if extension.eq_ignore_ascii_case("toml") => toml::from_str(&text)?,
        _ => return Err("polymer role profiles must use a .json or .toml suffix".into()),
    };
    Ok(PolymerRoleProfile {
        id: document.profile_id,
        rules: document
            .rules
            .into_iter()
            .map(|rule| PolymerRoleRule {
                atom_name: rule.atom_name,
                role: rule.role,
                component_id: rule.component_id,
                component_kind: rule.component_kind,
            })
            .collect(),
    })
}

/// Reads the caller's role profile, reporting a failure as a policy error.
pub(crate) fn load(profile: &Path) -> Result<PolymerRoleProfile, crate::exit::Exit> {
    read(profile).map_err(|error| {
        eprintln!("polymer role profile failed: {error}");
        crate::exit::Exit::Policy
    })
}

/// Applies a role profile and returns the role-annotated structure.
///
/// Backbone and nucleotide checks name atoms by role, never by atom name, so
/// the profile is the caller's explicit statement of which atom plays which
/// part under the given dictionary.
pub(crate) fn apply(
    structure: &molframe::Structure,
    profile: &PolymerRoleProfile,
    source: &Path,
    provider: &impl molframe::chemistry::ComponentProvider,
    context: crate::report::Context,
) -> Result<molframe_core::Structure, crate::exit::Exit> {
    match molframe::chemistry::apply_polymer_role_profile(structure.engine(), provider, profile) {
        Ok(report) => {
            for component in &report.unresolved_components {
                eprintln!(
                    "polymer role profile {}: unresolved CCD component {component}",
                    report.profile_id
                );
            }
            Ok(report.structure)
        }
        Err(finding) => {
            context.findings(&[finding], &source.display().to_string());
            Err(crate::exit::Exit::Policy)
        }
    }
}

#[cfg(test)]
#[path = "role_profile_tests.rs"]
mod tests;
