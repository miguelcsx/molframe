//! CCD-expected atom coverage under an explicit analysis policy.

use crate::ComponentProvider;
use pdbiox_core::contract::{AnalysisPolicy, Coverage, DictionaryVersion, HydrogenPolicy};
use pdbiox_core::{Code, Diagnostic, Element, Structure};
use std::collections::{BTreeMap, BTreeSet};

/// Coverage and traceability for one structure against component definitions.
#[derive(Clone, Debug)]
pub struct ComponentCoverage {
    /// Expected, used, missing and ambiguous atom counts.
    pub coverage: Coverage,
    /// Unknown components and policy-related findings.
    pub findings: Vec<Diagnostic>,
    /// Exact component dictionary version used.
    pub dictionary_version: DictionaryVersion,
}

/// Measures modelled atoms against the CCD-expected non-leaving atom set.
///
/// Alternate locations and hydrogen intent are resolved through `policy`.
/// Components absent from the provider are reported and excluded because their
/// intended count is unknowable; they are never treated as complete.
///
/// # Errors
///
/// Returns a provider diagnostic or an arithmetic-capacity diagnostic.
pub fn component_coverage(
    structure: &Structure,
    provider: &dyn ComponentProvider,
    policy: &AnalysisPolicy,
) -> Result<ComponentCoverage, Diagnostic> {
    let selected = structure.resolve_altlocs(policy);
    let mut findings = selected.warnings;
    let mut unresolved = BTreeSet::new();
    let mut coverage = Coverage::default();
    for residue in structure.data().residues() {
        let Some(component_id) = residue.name() else {
            continue;
        };
        let Some(component) = provider.get(component_id)? else {
            unresolved.insert(Box::<str>::from(component_id));
            continue;
        };
        let mut observed = BTreeMap::<&str, u32>::new();
        for atom in residue
            .atoms()
            .filter(|atom| selected.value.contains(atom.index().get()))
        {
            if atom.component_name() != Some(component_id) {
                continue;
            }
            let Some(name) = atom.name() else {
                continue;
            };
            let count = observed.entry(name).or_default();
            *count = count.checked_add(1).ok_or_else(capacity)?;
        }
        for expected in component.atoms.iter().filter(|atom| !atom.leaving) {
            let count = match observed.get(expected.name.as_ref()) {
                Some(count) => *count,
                None => 0,
            };
            if !hydrogen_is_intended(expected.element, count, policy.hydrogens) {
                continue;
            }
            coverage.intended = coverage.intended.checked_add(1).ok_or_else(capacity)?;
            match count {
                0 => coverage.missing = coverage.missing.checked_add(1).ok_or_else(capacity)?,
                1 => coverage.used = coverage.used.checked_add(1).ok_or_else(capacity)?,
                _ => {
                    coverage.ambiguous = coverage.ambiguous.checked_add(1).ok_or_else(capacity)?;
                }
            }
        }
    }
    findings.extend(
        unresolved
            .into_iter()
            .map(|component| Diagnostic::new(Code::W3201).with_context("component", component)),
    );
    Ok(ComponentCoverage {
        coverage,
        findings,
        dictionary_version: provider.version().clone(),
    })
}

fn hydrogen_is_intended(element: Element, observed: u32, policy: HydrogenPolicy) -> bool {
    if element != Element::HYDROGEN {
        return true;
    }
    if matches!(policy, HydrogenPolicy::IncludeInferred) {
        true
    } else if matches!(policy, HydrogenPolicy::ExplicitOnly) {
        observed > 0
    } else {
        false
    }
}

fn capacity() -> Diagnostic {
    Diagnostic::new(Code::E1901).with_context("limit", "component coverage count")
}

#[cfg(test)]
#[path = "coverage_tests.rs"]
mod tests;
