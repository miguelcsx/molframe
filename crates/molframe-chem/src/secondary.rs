//! Shared native DSSP 4 secondary-structure assignment.
//!
//! A sorted C-alpha cell grid bounds hydrogen-bond candidates to 9 angstroms.
//! For R residues and P candidate pairs, grid/pair sorting costs O(R log R +
//! P log P), then classification is local-density bounded. Storage is O(R + P).
//! No external executable is required. File precedence remains the caller's policy.

use crate::grid::{CellGrid, cell_for};
use molframe_core::hashing::IdentityHashSet;
use molframe_core::structure::Structure;
use molframe_core::{SecondaryAssignment, SecondarySource, SecondaryStructure as Ss};

#[path = "secondary/bridges.rs"]
mod bridges;
#[path = "secondary/classify.rs"]
mod classify;
#[path = "secondary/geometry.rs"]
mod geometry;
#[path = "secondary/options.rs"]
mod options;
use geometry::{bond_energy, squared_distance};
pub use options::{DsspBackbone, DsspOptions, InvalidDsspOptions};
type Bonds = IdentityHashSet<(usize, usize)>;
const CA_CUTOFF: f32 = 9.0;

/// Assigns native DSSP states from complete named heavy backbones.
///
/// Deposited records are not consulted. Missing heavy backbones remain unknown;
/// chains containing only C-alpha coordinates use a coarse trace fallback.
#[must_use]
pub fn assign_secondary_structure(structure: &Structure) -> Vec<SecondaryAssignment> {
    let mut backbones = vec![DsspBackbone::default(); structure.residue_count()];
    for model in structure.data().models() {
        for chain in model.chains() {
            for residue in chain.residues() {
                let slot = &mut backbones[residue.index().as_usize()];
                slot.model = model.index().get();
                slot.chain = chain.index().get();
                slot.proline = residue.name() == Some("PRO");
                for atom in residue.atoms() {
                    let target = match atom.name() {
                        Some("CA") => &mut slot.ca,
                        Some("N") => &mut slot.nitrogen,
                        Some("C") => &mut slot.carbon,
                        Some("O") => &mut slot.oxygen,
                        _ => continue,
                    };
                    if target.is_none() {
                        *target = atom.position();
                    }
                }
            }
        }
    }
    let mut states = classify_backbones(&backbones, &DsspOptions::default());
    let mut sources: Vec<_> = states
        .iter()
        .map(|s| {
            if *s == Ss::Unknown {
                SecondarySource::None
            } else {
                SecondarySource::Dssp
            }
        })
        .collect();
    for chain in structure.data().chains() {
        let range = structure.data().topology.chains.residues(chain.index());
        let Some(range) = range else {
            continue;
        };
        let first = range.start as usize;
        let last = range.end as usize;
        let trace = &backbones[first..last];
        if trace
            .iter()
            .any(|b| b.nitrogen.is_some() || b.carbon.is_some() || b.oxygen.is_some())
        {
            continue;
        }
        for i in first..last {
            if backbones[i].ca.is_some() {
                states[i] = Ss::Coil;
                sources[i] = SecondarySource::CaOnly;
            }
        }
        assign_ca_trace(trace, &mut states[first..last]);
    }
    states
        .into_iter()
        .zip(sources)
        .map(|(state, source)| SecondaryAssignment { state, source })
        .collect()
}

/// Classifies supplied backbones using the same kernel as automatic enrichment.
///
/// The input order is polymer order, with model/chain identities supplied by the
/// extractor. Missing roles remain `Unknown`, including missing carbonyl oxygen.
///
/// # Errors
///
/// Returns `InvalidDsspOptions` instead of inventing assignments for invalid parameters.
pub fn dssp_from_backbones(
    backbones: &[DsspBackbone],
    options: &DsspOptions,
) -> Result<Vec<Ss>, InvalidDsspOptions> {
    if !options.is_valid() {
        return Err(InvalidDsspOptions);
    }
    Ok(classify_backbones(backbones, options))
}

fn classify_backbones(backbones: &[DsspBackbone], options: &DsspOptions) -> Vec<Ss> {
    let mut states = backbones
        .iter()
        .map(|b| {
            if b.is_evaluable() {
                Ss::Coil
            } else {
                Ss::Unknown
            }
        })
        .collect::<Vec<_>>();
    let (pairs, bonds) = hydrogen_bonds(backbones, options);
    classify::classify(backbones, &bonds, &pairs, &mut states, options);
    for (backbone, state) in backbones.iter().zip(&mut states) {
        if !backbone.is_evaluable() {
            *state = Ss::Unknown;
        }
    }
    states
}

fn hydrogen_bonds(
    backbones: &[DsspBackbone],
    options: &DsspOptions,
) -> (Vec<(usize, usize)>, Bonds) {
    let mut entries = Vec::with_capacity(backbones.len());
    for (residue, backbone) in backbones.iter().enumerate() {
        if !backbone.is_evaluable() {
            continue;
        }
        let (Some(position), Ok(residue)) = (backbone.ca, u32::try_from(residue)) else {
            continue;
        };
        if let Some(cell) = cell_for(position, CA_CUTOFF) {
            entries.push((cell, residue));
        }
    }
    let grid = CellGrid::build(entries);
    let mut pairs = Vec::new();
    grid.for_each_cell(|own, neighbourhood| {
        for &left in &grid.items()[own] {
            for range in neighbourhood {
                for &right in &grid.items()[range.clone()] {
                    if right <= left {
                        continue;
                    }
                    let (i, j) = (left as usize, right as usize);
                    if backbones[i].model == backbones[j].model
                        && squared_distance(backbones[i].ca, backbones[j].ca)
                            <= CA_CUTOFF * CA_CUTOFF
                    {
                        pairs.push((i, j));
                    }
                }
            }
        }
    });
    pairs.sort_unstable();
    pairs.dedup();
    // Only the two strongest acceptors of each donor define DSSP TestBond.
    // Sorted residue pairs also make energy ties independent of grid traversal.
    let mut strongest = vec![[(usize::MAX, 0.0); 2]; backbones.len()];
    for &(i, j) in &pairs {
        if j - i < options.minimum_sequence_separation {
            continue;
        }
        let directions = [(j, i), (i, j)];
        for (acceptor, donor) in directions {
            if donor == j && j == i + 1 {
                continue;
            }
            let energy = bond_energy(backbones, acceptor, donor, options);
            let slots = &mut strongest[donor];
            if energy < slots[0].1 {
                slots[1] = slots[0];
                slots[0] = (acceptor, energy);
            } else if energy < slots[1].1 {
                slots[1] = (acceptor, energy);
            }
        }
    }
    let mut bonds = Bonds::default();
    for (donor, slots) in strongest.iter().enumerate() {
        for &(acceptor, energy) in slots {
            if acceptor != usize::MAX && energy < options.hydrogen_bond_energy {
                bonds.insert((acceptor, donor));
            }
        }
    }
    (pairs, bonds)
}

fn assign_ca_trace(backbones: &[DsspBackbone], states: &mut [Ss]) {
    for i in 0..backbones.len().saturating_sub(3) {
        let window = &backbones[i..i + 4];
        if window
            .windows(2)
            .any(|p| squared_distance(p[0].ca, p[1].ca) > 4.5 * 4.5)
        {
            continue;
        }
        let distance = squared_distance(window[0].ca, window[3].ca);
        if (3.08 * 3.08..=7.28 * 7.28).contains(&distance) {
            states[i..i + 4].fill(Ss::AlphaHelix);
        } else if (8.98 * 8.98..=11.82 * 11.82).contains(&distance) {
            for state in &mut states[i..i + 4] {
                if *state == Ss::Coil {
                    *state = Ss::Strand;
                }
            }
        }
    }
}

#[cfg(test)]
#[path = "secondary_tests.rs"]
mod tests;
