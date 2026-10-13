use crate::{AtomSite, ComponentRole, ComponentSpec, MappingPlan, Motif};
use molframe_chem::{Component, ComponentKind, ComponentProvider, equivalence_classes};
use molframe_core::contract::{AnalysisPolicy, EquivalencePolicy};
use molframe_core::execution::CancellationToken;
use molframe_core::index::{AtomIndex, ResidueIndex};
use molframe_core::structure::{ResidueRef, Structure};
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
    Chemistry(molframe_core::Diagnostic),
    /// The explicit enumeration bound was reached.
    #[error(
        "motif mapping exceeded the limit of {limit} mappings; raise the mapping limit, \
         or stream and count the mappings instead of materialising them"
    )]
    LimitExceeded {
        /// Caller-selected upper bound on the number of mappings.
        limit: usize,
    },
    /// The caller cancelled the enumeration.
    #[error("motif mapping was cancelled")]
    Cancelled,
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
    let plan = std::sync::Arc::new(MappingPlan::build(
        structure, motif, provider, policy, false,
    )?);
    let cancel = CancellationToken::new();
    let mut stream = plan.stream();
    let mut tuple = vec![0u32; plan.width()];
    let mut mappings = Vec::new();
    while stream.next_into(&mut tuple, &cancel)? {
        if mappings.len() == limit {
            return Err(MappingError::LimitExceeded { limit });
        }
        mappings.push(plan.mapping(&tuple));
    }
    Ok(MappingSet {
        ambiguous: mappings.len() > 1,
        mappings,
    })
}

#[derive(Clone, Debug)]
pub(crate) struct Candidate {
    pub(crate) residue: ResidueIndex,
    pub(crate) atoms: BTreeMap<AtomSite, Vec<AtomIndex>>,
}

pub(crate) fn component_candidates(
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
        .map(molframe_core::structure::AtomRef::index)
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

pub(crate) fn referenced_sites(motif: &Motif) -> BTreeSet<AtomSite> {
    motif
        .constraints()
        .iter()
        .flat_map(|named| named.constraint.sites().into_iter().cloned())
        .collect()
}

#[cfg(test)]
#[path = "mapping_tests.rs"]
mod tests;
