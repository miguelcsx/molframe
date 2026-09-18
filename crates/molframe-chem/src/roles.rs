//! Explicit, versioned assignment of semantic polymer atom roles.

use crate::{ComponentKind, ComponentProvider, PolymerAtomRole};
use molframe_core::contract::DictionaryVersion;
use molframe_core::{AnnotationColumn, AtomAnnotation, Code, Diagnostic, Presence, Structure};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

/// One data-driven component-atom role assignment.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PolymerRoleRule {
    /// Optional exact CCD component identifier constraint.
    pub component_id: Option<Box<str>>,
    /// Optional broad CCD component-kind constraint.
    pub component_kind: Option<ComponentKind>,
    /// CCD-local atom identifier.
    pub atom_name: Box<str>,
    /// Semantic role bits assigned when all constraints match.
    pub role: PolymerAtomRole,
}

/// Caller-selected, inspectable polymer atom-role profile.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PolymerRoleProfile {
    /// Stable profile identifier/version.
    pub id: Box<str>,
    /// Ordered role rules; matching roles are combined.
    pub rules: Arc<[PolymerRoleRule]>,
}

impl PolymerRoleProfile {
    /// Resolves combined roles for one CCD-local atom under this profile.
    #[must_use]
    pub fn role_for(&self, component: &crate::Component, atom_name: &str) -> PolymerAtomRole {
        self.rules
            .iter()
            .filter(|rule| rule_matches(rule, &component.id, component.kind, atom_name))
            .fold(PolymerAtomRole::UNKNOWN, |roles, rule| {
                roles.union(rule.role)
            })
    }
}

/// Result of applying one role profile to a structure.
#[derive(Clone, Debug)]
pub struct PolymerRoleReport {
    /// Snapshot containing the typed role annotation.
    pub structure: Structure,
    /// Components that the supplied CCD provider could not resolve.
    pub unresolved_components: Vec<Box<str>>,
    /// Exact dictionary version consulted.
    pub dictionary_version: DictionaryVersion,
    /// Exact role profile identifier.
    pub profile_id: Box<str>,
}

/// Applies an explicit role profile using CCD component identity and kind.
///
/// No standard atom names are built into this operation. Generic names such as
/// protein backbone atoms belong in a versioned profile supplied by the caller.
/// Rules may constrain an exact component, a component kind, or both.
///
/// # Errors
///
/// Returns a diagnostic for an empty/invalid profile or provider failure.
pub fn apply_polymer_role_profile(
    structure: &Structure,
    provider: &dyn ComponentProvider,
    profile: &PolymerRoleProfile,
) -> Result<PolymerRoleReport, Diagnostic> {
    validate_profile(profile)?;
    let mut data = structure.data().clone();
    let mut cache = BTreeMap::new();
    let mut unresolved = BTreeSet::new();
    let mut entries = Vec::with_capacity(structure.atom_count() as usize);
    for atom in structure.data().atoms() {
        let Some(component_id) = atom.component_name() else {
            entries.push((PolymerAtomRole::UNKNOWN.code(), Presence::Unknown));
            continue;
        };
        let component = match cache.get(component_id) {
            Some(component) => Some(Arc::clone(component)),
            None => provider.get(component_id)?,
        };
        let Some(component) = component else {
            unresolved.insert(Box::<str>::from(component_id));
            entries.push((PolymerAtomRole::UNKNOWN.code(), Presence::Unknown));
            continue;
        };
        cache
            .entry(Box::<str>::from(component_id))
            .or_insert_with(|| Arc::clone(&component));
        let Some(atom_name) = atom.name() else {
            entries.push((PolymerAtomRole::UNKNOWN.code(), Presence::Unknown));
            continue;
        };
        let role = profile.role_for(&component, atom_name);
        let presence = if role == PolymerAtomRole::UNKNOWN {
            Presence::Inapplicable
        } else {
            Presence::Present
        };
        entries.push((role.code(), presence));
    }
    let column = AnnotationColumn::from_entries(entries)
        .map_err(|error| Diagnostic::new(Code::E6009).with_context("cause", error.to_string()))?;
    let _ = data.annotations.insert(
        molframe_core::POLYMER_ATOM_ROLE_ANNOTATION,
        AtomAnnotation::Integer(column),
    );
    Ok(PolymerRoleReport {
        structure: Structure::new(data),
        unresolved_components: unresolved.into_iter().collect(),
        dictionary_version: provider.version().clone(),
        profile_id: profile.id.clone(),
    })
}

fn validate_profile(profile: &PolymerRoleProfile) -> Result<(), Diagnostic> {
    let invalid_rule = profile.rules.iter().any(|rule| {
        (rule.component_id.is_none() && rule.component_kind.is_none())
            || rule.atom_name.is_empty()
            || rule.role == PolymerAtomRole::UNKNOWN
    });
    if profile.id.is_empty() || profile.rules.is_empty() || invalid_rule {
        Err(Diagnostic::new(Code::E4003)
            .with_context("required", "non-empty versioned polymer role rules"))
    } else {
        Ok(())
    }
}

fn rule_matches(
    rule: &PolymerRoleRule,
    component_id: &str,
    component_kind: ComponentKind,
    atom_name: &str,
) -> bool {
    rule.atom_name.as_ref() == atom_name
        && rule
            .component_id
            .as_deref()
            .is_none_or(|required| required == component_id)
        && rule
            .component_kind
            .is_none_or(|required| required == component_kind)
}

#[cfg(test)]
#[path = "roles_tests.rs"]
mod tests;
