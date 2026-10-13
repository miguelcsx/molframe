//! Shared native DSSP 4 secondary-structure assignment.
//!
//! A sorted C-alpha cell grid bounds hydrogen-bond candidates to 9 angstroms.
//! For R residues and P candidate pairs, grid/pair sorting costs O(R log R +
//! P log P), then classification is local-density bounded. Storage is O(R + P).
//! No external executable is required. File precedence remains the caller's policy.

use crate::grid::{CellGrid, cell_for};
use molframe_core::ExecutionContext;
use molframe_core::hashing::IdentityHashSet;
use molframe_core::parallel::{BlockPlan, map_blocks_in};
use molframe_core::structure::ResidueRef;
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
    assign_secondary_structure_in(structure, &ExecutionContext::default())
}

/// [`assign_secondary_structure`] on the worker budget of `context`.
///
/// Backbone extraction, hydrogen-bond candidates, energies and the bridge
/// pre-filter run in fixed blocks merged in block order, so the result is
/// identical at every worker count.
#[must_use]
pub fn assign_secondary_structure_in(
    structure: &Structure,
    context: &ExecutionContext,
) -> Vec<SecondaryAssignment> {
    let backbones = gather_backbones(structure, context);
    let mut states = classify_backbones(&backbones, &DsspOptions::default(), context);
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
    Ok(classify_backbones(
        backbones,
        options,
        &ExecutionContext::default(),
    ))
}

fn classify_backbones(
    backbones: &[DsspBackbone],
    options: &DsspOptions,
    context: &ExecutionContext,
) -> Vec<Ss> {
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
    let (pairs, bonds) = hydrogen_bonds(backbones, options, context);
    classify::classify(backbones, &bonds, &pairs, &mut states, options, context);
    for (backbone, state) in backbones.iter().zip(&mut states) {
        if !backbone.is_evaluable() {
            *state = Ss::Unknown;
        }
    }
    states
}

const CHAINS_PER_BLOCK: usize = 4;
const RESIDUES_PER_BLOCK: usize = 8_192;
const PAIRS_PER_BLOCK: usize = 32_768;
const CELLS_PER_BLOCK: usize = 4096;

fn gather_backbones(structure: &Structure, context: &ExecutionContext) -> Vec<DsspBackbone> {
    let mut backbones = vec![DsspBackbone::default(); structure.residue_count()];
    let mut layout = Vec::new();
    for model in structure.data().models() {
        for chain in model.chains() {
            layout.push((model.index().get(), chain.index()));
        }
    }
    let plan = BlockPlan::new(layout.len(), CHAINS_PER_BLOCK);
    let run = |_: usize, range: std::ops::Range<usize>| {
        let mut out = Vec::new();
        for &(model, chain) in &layout[range] {
            let Some(chain) = structure.data().chain(chain) else {
                continue;
            };
            for residue in chain.residues() {
                out.push((
                    residue.index().as_usize(),
                    residue_backbone(model, chain.index().get(), residue),
                ));
            }
        }
        out
    };
    let blocks = match map_blocks_in(plan, context, run) {
        Ok(blocks) => blocks,
        // A worker failure falls back to the serial walk rather than guessing.
        Err(_) => vec![run(0, 0..layout.len())],
    };
    for (index, backbone) in blocks.into_iter().flatten() {
        backbones[index] = backbone;
    }
    backbones
}

fn residue_backbone(model: u32, chain: u32, residue: ResidueRef<'_>) -> DsspBackbone {
    let mut slot = DsspBackbone {
        model,
        chain,
        proline: residue.name() == Some("PRO"),
        ..DsspBackbone::default()
    };
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
    slot
}

fn hydrogen_bonds(
    backbones: &[DsspBackbone],
    options: &DsspOptions,
    context: &ExecutionContext,
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
    let search = |cells: std::ops::Range<usize>| {
        let mut found = Vec::new();
        grid.visit_cells(cells, |own, neighbourhood| {
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
                            found.push((i, j));
                        }
                    }
                }
            }
        });
        found
    };
    let cells = BlockPlan::new(grid.cell_count(), CELLS_PER_BLOCK);
    let blocks = match map_blocks_in(cells, context, |_, range| search(range)) {
        Ok(blocks) => blocks,
        Err(_) => vec![search(0..grid.cell_count())],
    };
    let pairs = sorted_pairs(blocks, backbones.len(), context);
    // Only the two strongest acceptors of each donor define DSSP TestBond.
    // Sorted residue pairs also make energy ties independent of grid traversal.
    // The energies are pure, so they are computed in blocks; the strongest-two
    // update then folds them serially in sorted pair order.
    let energy_of = |range: std::ops::Range<usize>| {
        pairs[range]
            .iter()
            .map(|&(i, j)| {
                if j - i < options.minimum_sequence_separation {
                    return [f64::INFINITY; 2];
                }
                [
                    // Acceptor j, donor i.
                    bond_energy(backbones, j, i, options),
                    // Acceptor i, donor j; adjacent residues have no such bond.
                    if j == i + 1 {
                        f64::INFINITY
                    } else {
                        bond_energy(backbones, i, j, options)
                    },
                ]
            })
            .collect::<Vec<_>>()
    };
    let plan = BlockPlan::new(pairs.len(), PAIRS_PER_BLOCK);
    let energies = match map_blocks_in(plan, context, |_, range| energy_of(range)) {
        Ok(energies) => energies,
        Err(_) => vec![energy_of(0..pairs.len())],
    };
    let mut strongest = vec![[(usize::MAX, 0.0_f64); 2]; backbones.len()];
    let mut position = 0;
    for block in &energies {
        for &[to_i, to_j] in block {
            let (i, j) = pairs[position];
            position += 1;
            for (acceptor, donor, energy) in [(j, i, to_i), (i, j, to_j)] {
                let slots = &mut strongest[donor];
                if energy < slots[0].1 {
                    slots[1] = slots[0];
                    slots[0] = (acceptor, energy);
                } else if energy < slots[1].1 {
                    slots[1] = (acceptor, energy);
                }
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

/// Merges per-block candidate pairs into one list sorted by `(i, j)` without
/// duplicates: a counting pass by `i`, then each `i` row sorted in blocks.
fn sorted_pairs(
    blocks: Vec<Vec<(usize, usize)>>,
    residues: usize,
    context: &ExecutionContext,
) -> Vec<(usize, usize)> {
    let mut starts = vec![0_usize; residues + 1];
    for &(i, _) in blocks.iter().flatten() {
        starts[i + 1] += 1;
    }
    for i in 0..residues {
        starts[i + 1] += starts[i];
    }
    let mut fill = starts.clone();
    let mut partners = vec![0_usize; starts[residues]];
    for &(i, j) in blocks.iter().flatten() {
        partners[fill[i]] = j;
        fill[i] += 1;
    }
    drop(blocks);
    let rows = BlockPlan::new(residues, RESIDUES_PER_BLOCK);
    let sort_rows = |range: std::ops::Range<usize>| {
        let mut out = Vec::with_capacity(starts[range.end] - starts[range.start]);
        for i in range {
            let mut row = partners[starts[i]..starts[i + 1]].to_vec();
            row.sort_unstable();
            row.dedup();
            out.extend(row.into_iter().map(|j| (i, j)));
        }
        out
    };
    let sorted = match map_blocks_in(rows, context, |_, range| sort_rows(range)) {
        Ok(sorted) => sorted,
        Err(_) => vec![sort_rows(0..residues)],
    };
    sorted.into_iter().flatten().collect()
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
