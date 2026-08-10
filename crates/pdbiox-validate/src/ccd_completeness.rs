//! Per-residue missing-atom reports against explicit CCD chemistry.

use pdbiox_chem::ComponentProvider;
use pdbiox_core::contract::{AnalysisPolicy, HydrogenPolicy, Namespace};
use pdbiox_core::index::ResidueIndex;
use pdbiox_core::structure::{AtomRef, ResidueRef, Structure};
use pdbiox_core::{Code, Diagnostic, Element};
use std::collections::BTreeMap;

/// CCD completeness for one topology residue.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ResidueAtomCompleteness {
    /// Topology residue index.
    pub residue: ResidueIndex,
    /// CCD component identifier.
    pub component: String,
    /// Expected non-leaving atoms under hydrogen policy.
    pub intended: usize,
    /// Expected atoms observed exactly once.
    pub assessed: usize,
    /// Missing CCD atom identifiers.
    pub missing: Vec<String>,
    /// CCD atom identifiers observed more than once after altloc resolution.
    pub ambiguous: Vec<String>,
}

/// Exact coverage and residue-local missing names.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CcdCompletenessReport {
    /// Expected atom count.
    pub intended: usize,
    /// Exactly-once observed expected atom count.
    pub assessed: usize,
    /// Multiply observed expected atom count.
    pub ambiguous: usize,
    /// Residue reports in topology order.
    pub residues: Vec<ResidueAtomCompleteness>,
}

/// Reports missing CCD atoms for every residue without compatibility tables.
///
/// Missing CCD components and identifiers are errors: this validator never
/// substitutes local residue-name chemistry.
///
/// # Errors
///
/// Returns provider, namespace, identifier or coverage diagnostics.
pub fn ccd_missing_atoms(
    structure: &Structure,
    provider: &dyn ComponentProvider,
    policy: &AnalysisPolicy,
) -> Result<CcdCompletenessReport, Diagnostic> {
    let selected = structure.resolve_altlocs(policy).value;
    let mut report = CcdCompletenessReport {
        intended: 0,
        assessed: 0,
        ambiguous: 0,
        residues: Vec::new(),
    };
    for residue in structure.data().residues() {
        let component_id =
            residue_name(residue, policy.identifiers)?.ok_or_else(identifier_error)?;
        let component = provider
            .get(component_id)?
            .ok_or_else(|| missing_component(component_id))?;
        let observed = observed_atoms(residue, &selected, policy.identifiers)?;
        let mut local = ResidueAtomCompleteness {
            residue: residue.index(),
            component: component_id.to_string(),
            intended: 0,
            assessed: 0,
            missing: Vec::new(),
            ambiguous: Vec::new(),
        };
        for expected in component.atoms.iter().filter(|atom| !atom.leaving) {
            let count = match observed.get(expected.name.as_ref()).copied() {
                Some(count) => count,
                None => 0,
            };
            if !hydrogen_intended(expected.element, count, policy.hydrogens) {
                continue;
            }
            local.intended += 1;
            match count {
                0 => local.missing.push(expected.name.to_string()),
                1 => local.assessed += 1,
                _ => local.ambiguous.push(expected.name.to_string()),
            }
        }
        report.intended += local.intended;
        report.assessed += local.assessed;
        report.ambiguous += local.ambiguous.len();
        report.residues.push(local);
    }
    Ok(report)
}

fn observed_atoms<'a>(
    residue: ResidueRef<'a>,
    selected: &pdbiox_core::AtomSelection,
    namespace: Namespace,
) -> Result<BTreeMap<&'a str, usize>, Diagnostic> {
    let mut observed = BTreeMap::new();
    for atom in residue
        .atoms()
        .filter(|atom| selected.contains(atom.index().get()))
    {
        let name = atom_name(atom, namespace)?.ok_or_else(identifier_error)?;
        *observed.entry(name).or_default() += 1;
    }
    Ok(observed)
}

fn residue_name(residue: ResidueRef<'_>, namespace: Namespace) -> Result<Option<&str>, Diagnostic> {
    match namespace {
        Namespace::Label => Ok(residue.name()),
        Namespace::Auth => Ok(residue.auth_name()),
        _ => Err(namespace_error()),
    }
}

fn atom_name(atom: AtomRef<'_>, namespace: Namespace) -> Result<Option<&str>, Diagnostic> {
    match namespace {
        Namespace::Label => Ok(atom.name()),
        Namespace::Auth => Ok(atom.auth_name()),
        _ => Err(namespace_error()),
    }
}

fn hydrogen_intended(element: Element, observed: usize, policy: HydrogenPolicy) -> bool {
    element != Element::HYDROGEN
        || matches!(policy, HydrogenPolicy::IncludeInferred)
        || matches!(policy, HydrogenPolicy::ExplicitOnly) && observed > 0
}

fn namespace_error() -> Diagnostic {
    Diagnostic::new(Code::E4003).with_context("namespace", "label or auth required")
}

fn identifier_error() -> Diagnostic {
    Diagnostic::new(Code::E4003).with_context("required", "component and atom identifiers")
}

fn missing_component(component: &str) -> Diagnostic {
    Diagnostic::new(Code::E4003)
        .with_context("required", "CCD component")
        .with_context("component", component)
}
