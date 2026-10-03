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

pub(super) fn read(path: &Path) -> Result<PolymerRoleProfile, Box<dyn std::error::Error>> {
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

#[cfg(test)]
#[path = "role_profile_tests.rs"]
mod tests;
