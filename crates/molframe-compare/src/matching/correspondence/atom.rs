//! Atom correspondence inside residue pairs.
//!
//! Atoms match by name. Where the dictionary defines the residue, chemically
//! equivalent atoms (the two carbons of a phenyl ring, the oxygens of a
//! carboxylate) may swap names between structures; the swap that minimises the
//! squared distance between mapped atoms is taken, and atoms a structure lacks
//! (typically hydrogens) do not vote. Distances are measured in the coordinate
//! frame the inputs are already in, so structures compared this way should
//! share a frame.

use super::residue::ResiduePair;
use crate::mapped::MappedCompareError;
use crate::workflow::{PointMapping, PointMatch};
use crate::{component_fragment, equivalent_atom_mappings};
use molframe_chem::{Component, ComponentProvider};
use molframe_core::contract::Status;
use molframe_core::index::AtomIndex;
use molframe_core::structure::{ResidueRef, Structure};
use std::collections::BTreeMap;

/// The atoms that correspond across two structures.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AtomCorrespondence {
    /// Reference atom to target atom, in reference residue order.
    pub mapping: PointMapping,
    /// Residue pairs whose equivalent atoms were swapped to fit better.
    pub swapped_residues: usize,
    /// Ambiguous when multiple minimum-distance CCD-equivalent mappings tie.
    pub status: Status,
    /// Equal-cost atom mappings, grouped by residue pair; each list includes
    /// the selected mapping and every tied alternative.
    pub alternatives: Vec<AtomMappingAlternative>,
}

/// Minimum-distance atom mappings for one aligned residue pair.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AtomMappingAlternative {
    /// The aligned residue indices in reference and target order.
    pub residues: ResiduePair,
    /// Tied mappings, expressed as global atom indices.
    pub mappings: Vec<Vec<PointMatch>>,
}

type AtomPairs = Vec<(AtomIndex, AtomIndex)>;

/// Matches atoms by name within each residue pair, resolving equivalences.
///
/// `automorphism_limit` bounds the equivalences enumerated per component; a
/// component with more is refused rather than truncated.
///
/// # Errors
///
/// Returns a provider failure, an exceeded automorphism bound, or an invalid
/// mapping when no atom corresponds.
pub fn map_atoms(
    reference: &Structure,
    target: &Structure,
    residues: &[ResiduePair],
    provider: &dyn ComponentProvider,
    automorphism_limit: usize,
) -> Result<AtomCorrespondence, MappedCompareError> {
    let mut matches = Vec::new();
    let mut swapped_residues = 0;
    let mut alternatives = Vec::new();
    for &(first, second) in residues {
        let (Some(first), Some(second)) = (
            reference.data().residue(first),
            target.data().residue(second),
        ) else {
            continue;
        };
        let component = match first.name() {
            Some(name) if first.name() == second.name() => {
                provider.get(name).map_err(MappedCompareError::Mapping)?
            }
            _ => None,
        };
        let pairs = match component {
            Some(component) => {
                let (pairs, swapped, tied) =
                    equivalent_pairs(first, second, &component, automorphism_limit)?;
                swapped_residues += usize::from(swapped);
                if tied.len() > 1 {
                    alternatives.push(AtomMappingAlternative {
                        residues: (first.index(), second.index()),
                        mappings: tied
                            .into_iter()
                            .map(|pairs| {
                                pairs
                                    .into_iter()
                                    .map(|(reference, model)| PointMatch {
                                        reference: reference.as_usize(),
                                        model: model.as_usize(),
                                    })
                                    .collect()
                            })
                            .collect(),
                    });
                }
                pairs
            }
            None => named_pairs(first, second),
        };
        matches.extend(pairs.into_iter().map(|(reference, target)| PointMatch {
            reference: reference.as_usize(),
            model: target.as_usize(),
        }));
    }
    let mapping = PointMapping::new(
        matches,
        reference.atom_count() as usize,
        target.atom_count() as usize,
    )
    .map_err(MappedCompareError::Compare)?;
    Ok(AtomCorrespondence {
        mapping,
        swapped_residues,
        status: if alternatives.is_empty() {
            Status::Complete
        } else {
            Status::Ambiguous
        },
        alternatives,
    })
}

/// First atom of each name in a residue.
fn atoms_by_name(residue: ResidueRef<'_>) -> BTreeMap<&str, AtomIndex> {
    let mut atoms = BTreeMap::new();
    for atom in residue.atoms() {
        if let Some(name) = atom.name() {
            atoms.entry(name).or_insert(atom.index());
        }
    }
    atoms
}

fn named_pairs(first: ResidueRef<'_>, second: ResidueRef<'_>) -> AtomPairs {
    let sources = atoms_by_name(first);
    let targets = atoms_by_name(second);
    // File order, so the mapping is stable.
    first
        .atoms()
        .filter_map(|atom| {
            let name = atom.name()?;
            (sources.get(name) == Some(&atom.index()))
                .then(|| targets.get(name).map(|target| (atom.index(), *target)))
                .flatten()
        })
        .collect()
}

fn position_of(residue: ResidueRef<'_>, atom: AtomIndex) -> Option<[f32; 3]> {
    residue
        .atoms()
        .find(|candidate| candidate.index() == atom)
        .and_then(molframe_core::structure::AtomRef::position)
}

/// Name-matched atoms, with equivalent atoms permuted to the nearest fit.
fn equivalent_pairs(
    first: ResidueRef<'_>,
    second: ResidueRef<'_>,
    component: &Component,
    limit: usize,
) -> Result<(AtomPairs, bool, Vec<AtomPairs>), MappedCompareError> {
    let sources = atoms_by_name(first);
    let targets = atoms_by_name(second);
    let component = &component_fragment(component, |name| {
        sources.contains_key(name) && targets.contains_key(name)
    });
    let cost = |mapping: &[u32]| -> f64 {
        let mut total = 0.0;
        for (index, atom) in component.atoms.iter().enumerate() {
            let Some(source) = sources.get(&*atom.name) else {
                continue;
            };
            let model_name = mapping
                .get(index)
                .and_then(|slot| component.atoms.get(usize::try_from(*slot).ok()?))
                .map(|target| &*target.name);
            let Some(target) = model_name.and_then(|name| targets.get(name)) else {
                return f64::INFINITY;
            };
            let (Some(from), Some(to)) =
                (position_of(first, *source), position_of(second, *target))
            else {
                return f64::INFINITY;
            };
            total += molframe_geom::distance_squared(from, to);
        }
        total
    };
    let mut best: Option<(f64, Box<[u32]>)> = None;
    let mut tied = Vec::new();
    for candidate in
        equivalent_atom_mappings(component, limit).map_err(MappedCompareError::Compare)?
    {
        let total = cost(&candidate.reference_to_model);
        let mapping = candidate.reference_to_model;
        match best.as_ref() {
            None => {
                best = Some((total, mapping.clone()));
                tied.push(mapping);
            }
            Some((current, _)) if strictly_better(total, *current) => {
                best = Some((total, mapping.clone()));
                tied.clear();
                tied.push(mapping);
            }
            Some((current, _)) if tied_equal(total, *current) => tied.push(mapping),
            _ => {}
        }
    }
    let Some((_, best_mapping)) = best else {
        return Ok((named_pairs(first, second), false, Vec::new()));
    };
    let (pairs, swapped) = pairs_for_mapping(component, &sources, &targets, &best_mapping);
    let alternatives = if tied.len() > 1 {
        tied.into_iter()
            .map(|mapping| pairs_for_mapping(component, &sources, &targets, &mapping).0)
            .collect()
    } else {
        Vec::new()
    };
    Ok((pairs, swapped, alternatives))
}

fn pairs_for_mapping(
    component: &Component,
    sources: &BTreeMap<&str, AtomIndex>,
    targets: &BTreeMap<&str, AtomIndex>,
    mapping: &[u32],
) -> (AtomPairs, bool) {
    let mut pairs = Vec::new();
    let mut swapped = false;
    for (index, atom) in component.atoms.iter().enumerate() {
        let target_name = mapping
            .get(index)
            .and_then(|slot| component.atoms.get(usize::try_from(*slot).ok()?))
            .map(|target| &*target.name);
        if let (Some(source), Some(target)) = (
            sources.get(&*atom.name),
            target_name.and_then(|name| targets.get(name)),
        ) {
            swapped |= target_name != Some(&*atom.name);
            pairs.push((*source, *target));
        }
    }
    for (name, source) in sources {
        if component.atom(name).is_none()
            && let Some(target) = targets.get(name)
        {
            pairs.push((*source, *target));
        }
    }
    (pairs, swapped)
}

fn tied_equal(left: f64, right: f64) -> bool {
    if !left.is_finite() || !right.is_finite() {
        return false;
    }
    let tolerance = molframe_core::contract::Tolerance::default();
    (left - right).abs() <= tolerance.absolute + tolerance.relative * left.abs().max(right.abs())
}

fn strictly_better(candidate: f64, current: f64) -> bool {
    candidate < current && !tied_equal(candidate, current)
}
