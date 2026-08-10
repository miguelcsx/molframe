use crate::{AtomSite, ComponentRole, ComponentSpec, Motif};
use pdbiox_chem::{Component, ComponentKind, ComponentProvider, equivalence_classes};
use pdbiox_core::contract::{AnalysisPolicy, EquivalencePolicy};
use pdbiox_core::index::{AtomIndex, ResidueIndex};
use pdbiox_core::structure::{ResidueRef, Structure};
use std::collections::{BTreeMap, BTreeSet};

use crate::numeric::usize_to_u32;

/// One complete component and atom correspondence.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MappedMotif {
    /// Motif component to structure residue.
    pub components: BTreeMap<Box<str>, ResidueIndex>,
    /// Every referenced motif atom and its chemically valid alternatives.
    pub atoms: BTreeMap<AtomSite, Vec<AtomIndex>>,
}

/// All defensible mappings in deterministic order.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MappingSet {
    /// Complete mappings.
    pub mappings: Vec<MappedMotif>,
    /// More than one complete mapping exists.
    pub ambiguous: bool,
}

/// Failure while enumerating motif mappings.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum MappingError {
    /// The CCD provider failed.
    #[error("component chemistry could not be read: {0}")]
    Chemistry(pdbiox_core::Diagnostic),
    /// The explicit enumeration bound was reached.
    #[error("motif mapping exceeded the limit of {limit} alternatives")]
    LimitExceeded {
        /// Caller-selected upper bound.
        limit: usize,
    },
}

/// Enumerates complete chemistry-aware motif mappings.
///
/// Component candidates are ordered by residue index. Under
/// [`EquivalencePolicy::Ccd`], atom alternatives come from exact CCD graph
/// automorphism orbits; explicitly declared equivalence groups are always
/// honoured. Residues cannot satisfy two motif components in one mapping.
///
/// # Errors
///
/// Returns [`MappingError::Chemistry`] when the provider fails and
/// [`MappingError::LimitExceeded`] before producing a partial enumeration.
pub fn map_motif(
    structure: &Structure,
    motif: &Motif,
    provider: Option<&dyn ComponentProvider>,
    policy: &AnalysisPolicy,
    limit: usize,
) -> Result<MappingSet, MappingError> {
    let sites = referenced_sites(motif);
    let mut candidates = BTreeMap::new();
    for (name, spec) in motif.components() {
        candidates.insert(
            name.clone(),
            component_candidates(structure, name, spec, provider, policy, &sites)?,
        );
    }
    let ordered: Vec<_> = motif.components().keys().cloned().collect();
    let mut state = EnumerationState::default();
    enumerate(&ordered, &candidates, 0, &mut state, limit)?;
    Ok(MappingSet {
        ambiguous: state.output.len() > 1,
        mappings: state.output,
    })
}

#[derive(Clone)]
struct Candidate {
    residue: ResidueIndex,
    atoms: BTreeMap<AtomSite, Vec<AtomIndex>>,
}

fn component_candidates(
    structure: &Structure,
    motif_name: &str,
    spec: &ComponentSpec,
    provider: Option<&dyn ComponentProvider>,
    policy: &AnalysisPolicy,
    sites: &BTreeSet<AtomSite>,
) -> Result<Vec<Candidate>, MappingError> {
    let mut candidates = Vec::new();
    for residue in structure.data().residues() {
        let Some(component_id) = residue.name() else {
            continue;
        };
        if !spec.component_ids.is_empty() && !spec.component_ids.contains(component_id) {
            continue;
        }
        let chemistry = match provider {
            Some(provider) => provider
                .get(component_id)
                .map_err(MappingError::Chemistry)?,
            None => None,
        };
        if !role_matches(
            spec.role,
            chemistry.as_deref(),
            !spec.component_ids.is_empty(),
        ) {
            continue;
        }
        let mut atoms = BTreeMap::new();
        let relevant = sites
            .iter()
            .filter(|site| site.component.as_ref() == motif_name);
        let mut complete = true;
        for site in relevant {
            let alternatives = atom_alternatives(residue, site, spec, chemistry.as_deref(), policy);
            if alternatives.is_empty() {
                complete = false;
                break;
            }
            atoms.insert(site.clone(), alternatives);
        }
        if complete
            && spec.required_atoms.iter().all(|required| {
                !atom_alternatives(
                    residue,
                    &AtomSite::new("", required.clone()),
                    spec,
                    chemistry.as_deref(),
                    policy,
                )
                .is_empty()
            })
        {
            candidates.push(Candidate {
                residue: residue.index(),
                atoms,
            });
        }
    }
    Ok(candidates)
}

fn role_matches(
    role: ComponentRole,
    chemistry: Option<&Component>,
    explicitly_classified: bool,
) -> bool {
    let Some(component) = chemistry else {
        return explicitly_classified;
    };
    match role {
        ComponentRole::Residue => matches!(
            component.kind,
            ComponentKind::AminoAcid | ComponentKind::Nucleotide | ComponentKind::Saccharide
        ),
        ComponentRole::Ligand | ComponentRole::Cofactor => matches!(
            component.kind,
            ComponentKind::NonPolymer | ComponentKind::Lipid
        ),
        ComponentRole::Metal => component.kind == ComponentKind::Ion,
    }
}

fn atom_alternatives(
    residue: ResidueRef<'_>,
    site: &AtomSite,
    spec: &ComponentSpec,
    chemistry: Option<&Component>,
    policy: &AnalysisPolicy,
) -> Vec<AtomIndex> {
    let mut names = BTreeSet::from([site.atom.clone()]);
    for group in &spec.equivalent_atoms {
        if group.contains(site.atom.as_ref()) {
            names.extend(group.iter().cloned());
        }
    }
    if policy.atom_equivalence == EquivalencePolicy::Ccd
        && let Some(component) = chemistry
    {
        expand_ccd_names(&mut names, component, &site.atom);
    }
    residue
        .atoms()
        .filter(|atom| atom.name().is_some_and(|name| names.contains(name)))
        .map(pdbiox_core::structure::AtomRef::index)
        .collect()
}

fn expand_ccd_names(names: &mut BTreeSet<Box<str>>, component: &Component, atom: &str) {
    let Some(position) = component
        .atoms
        .iter()
        .position(|candidate| candidate.name.as_ref() == atom)
    else {
        return;
    };
    let classes = equivalence_classes(component);
    let Some(class) = classes.class_of(usize_to_u32(position)) else {
        return;
    };
    for (index, candidate) in component.atoms.iter().enumerate() {
        if classes.class_of(usize_to_u32(index)) == Some(class) {
            names.insert(candidate.name.clone());
        }
    }
}

fn referenced_sites(motif: &Motif) -> BTreeSet<AtomSite> {
    motif
        .constraints()
        .iter()
        .flat_map(|named| named.constraint.sites().into_iter().cloned())
        .collect()
}

fn enumerate(
    ordered: &[Box<str>],
    candidates: &BTreeMap<Box<str>, Vec<Candidate>>,
    depth: usize,
    state: &mut EnumerationState,
    limit: usize,
) -> Result<(), MappingError> {
    if depth == ordered.len() {
        if state.output.len() == limit {
            return Err(MappingError::LimitExceeded { limit });
        }
        state.output.push(MappedMotif {
            components: state.components.clone(),
            atoms: state.atoms.clone(),
        });
        return Ok(());
    }
    let name = &ordered[depth];
    let Some(options) = candidates.get(name) else {
        return Ok(());
    };
    for candidate in options {
        if !state.used.insert(candidate.residue) {
            continue;
        }
        state.components.insert(name.clone(), candidate.residue);
        let old_atoms = state.atoms.clone();
        state.atoms.extend(candidate.atoms.clone());
        enumerate(ordered, candidates, depth + 1, state, limit)?;
        state.atoms = old_atoms;
        state.components.remove(name);
        state.used.remove(&candidate.residue);
    }
    Ok(())
}

#[derive(Default)]
struct EnumerationState {
    used: BTreeSet<ResidueIndex>,
    components: BTreeMap<Box<str>, ResidueIndex>,
    atoms: BTreeMap<AtomSite, Vec<AtomIndex>>,
    output: Vec<MappedMotif>,
}

#[cfg(test)]
#[path = "mapping_tests.rs"]
mod tests;
