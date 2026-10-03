//! Lightweight automatic secondary-structure assignment.
//!
//! The Kabsch–Sander definition of DSSP: a backbone hydrogen bond is an
//! electrostatic energy below −0.5 kcal/mol between a carbonyl and an amide,
//! helices are runs of two overlapping n-turns, and strands are residues in a
//! ladder of two or more consecutive bridges. A sorted cell grid limits the
//! candidate pairs to Cα pairs within 9 Å, so the pass is O(R·k) for R residues
//! and bounded local density.

use crate::grid::{CellGrid, cell_for};
use molframe_core::hashing::{IdentityHashMap, IdentityHashSet};
use molframe_core::structure::Structure;
use molframe_core::{SecondaryAssignment, SecondarySource, SecondaryStructure};

const CA_CUTOFF: f32 = 9.0;
/// The Kabsch–Sander bond threshold, in kcal/mol.
const HBOND_CUTOFF: f64 = -0.5;
/// The longest C–N distance that still counts as a peptide bond.
const PEPTIDE_BOND_LIMIT: f32 = 2.5;
/// Where the amide hydrogen sits along the previous carbonyl's C←O axis.
const NH_LENGTH: f32 = 1.0;
/// 0.084 · 332: the partial charges of the two dipoles times the unit factor.
const ELECTROSTATIC_FACTOR: f64 = 27.888;
/// DSSP's bend threshold on the change of Cα direction, in degrees.
const BEND_KAPPA_DEGREES: f64 = 70.0;

#[derive(Clone, Copy, Debug, Default)]
struct Backbone {
    chain: u32,
    proline: bool,
    ca: Option<[f32; 3]>,
    carbon: Option<[f32; 3]>,
    oxygen: Option<[f32; 3]>,
    nitrogen: Option<[f32; 3]>,
}

impl Backbone {
    const fn is_evaluable(&self) -> bool {
        self.ca.is_some() && self.nitrogen.is_some() && self.carbon.is_some()
    }
}

/// One carbonyl-to-amide hydrogen bond list, as `(acceptor, donor)` residues.
type Bonds = IdentityHashSet<(usize, usize)>;

/// Assigns every residue with a grid-bounded DSSP pass.
///
/// File records are intentionally not consulted here; the caller merges this
/// assignment with the file's by [`SecondarySource::rank`]. Residues with full
/// backbones get [`SecondarySource::Dssp`]; a Cα-only trace falls back to the
/// Zhang–Skolnick helix rule with [`SecondarySource::CaOnly`].
#[must_use]
pub fn assign_secondary_structure(structure: &Structure) -> Vec<SecondaryAssignment> {
    let mut backbones = vec![Backbone::default(); structure.residue_count()];
    for chain in structure.data().chains() {
        for residue in chain.residues() {
            let Some(slot) = backbones.get_mut(residue.index().as_usize()) else {
                continue;
            };
            slot.chain = chain.index().get();
            slot.proline = residue.name() == Some("PRO");
            // An alternate location repeats the atom name; the first
            // conformer is the one DSSP reads.
            for atom in residue.atoms() {
                let position = atom.position();
                let target = match atom.name() {
                    Some("CA") => &mut slot.ca,
                    Some("N") => &mut slot.nitrogen,
                    Some("C") => &mut slot.carbon,
                    Some("O") => &mut slot.oxygen,
                    _ => continue,
                };
                if target.is_none() {
                    *target = position;
                }
            }
        }
    }

    let mut states = backbones
        .iter()
        .map(|backbone| {
            if backbone.is_evaluable() {
                SecondaryStructure::Coil
            } else {
                SecondaryStructure::Unknown
            }
        })
        .collect::<Vec<_>>();
    let source = if backbones.iter().any(Backbone::is_evaluable) {
        let (pairs, bonds) = hydrogen_bonds(&backbones);
        classify(&backbones, &bonds, &pairs, &mut states);
        SecondarySource::Dssp
    } else {
        for (backbone, state) in backbones.iter().zip(&mut states) {
            if backbone.ca.is_some() {
                *state = SecondaryStructure::Coil;
            }
        }
        assign_zhang_skolnick(&backbones, &mut states);
        SecondarySource::CaOnly
    };
    states
        .into_iter()
        .map(|state| SecondaryAssignment {
            state,
            source: match state {
                SecondaryStructure::Unknown => SecondarySource::None,
                _ => source,
            },
        })
        .collect()
}

/// Every hydrogen bond between residues of one chain, and the residue pairs
/// close enough to have been considered.
fn hydrogen_bonds(backbones: &[Backbone]) -> (Vec<(usize, usize)>, Bonds) {
    let mut entries = Vec::new();
    for (residue, backbone) in backbones.iter().enumerate() {
        let (Some(position), Ok(residue)) = (backbone.ca, u32::try_from(residue)) else {
            continue;
        };
        let Some(cell) = cell_for(position, CA_CUTOFF) else {
            continue;
        };
        entries.push((cell, residue));
    }
    let grid = CellGrid::build(entries);
    let mut pairs = Vec::new();
    let mut bonds = Bonds::default();
    grid.for_each_cell(|own, neighbourhood| {
        for &left in &grid.items()[own] {
            for range in neighbourhood {
                for &right in &grid.items()[range.clone()] {
                    if right <= left {
                        continue;
                    }
                    let (left, right) = (left as usize, right as usize);
                    let (Some(first), Some(second)) = (backbones.get(left), backbones.get(right))
                    else {
                        continue;
                    };
                    if first.chain != second.chain
                        || squared_distance(first.ca, second.ca) > CA_CUTOFF * CA_CUTOFF
                    {
                        continue;
                    }
                    pairs.push((left, right));
                    if right - left < 3 {
                        continue;
                    }
                    if bond_energy(backbones, left, right) < HBOND_CUTOFF {
                        bonds.insert((left, right));
                    }
                    if bond_energy(backbones, right, left) < HBOND_CUTOFF {
                        bonds.insert((right, left));
                    }
                }
            }
        }
    });
    (pairs, bonds)
}

/// The Kabsch–Sander energy of the carbonyl of `acceptor` against the amide of
/// `donor`, in kcal/mol; infinite where the geometry is missing.
fn bond_energy(backbones: &[Backbone], acceptor: usize, donor: usize) -> f64 {
    let (Some(carbonyl), Some(amide)) = (backbones.get(acceptor), backbones.get(donor)) else {
        return f64::INFINITY;
    };
    let (Some(carbon), Some(oxygen), Some(nitrogen), Some(hydrogen)) = (
        carbonyl.carbon,
        carbonyl.oxygen,
        amide.nitrogen,
        amide_hydrogen(backbones, donor),
    ) else {
        return f64::INFINITY;
    };
    let on = distance(oxygen, nitrogen);
    let ch = distance(carbon, hydrogen);
    let oh = distance(oxygen, hydrogen);
    let cn = distance(carbon, nitrogen);
    if on <= f32::EPSILON || ch <= f32::EPSILON || oh <= f32::EPSILON || cn <= f32::EPSILON {
        return f64::INFINITY;
    }
    ELECTROSTATIC_FACTOR
        * (1.0 / f64::from(on) + 1.0 / f64::from(ch) - 1.0 / f64::from(oh) - 1.0 / f64::from(cn))
}

/// The amide hydrogen of `residue`, placed antiparallel to the C=O of the
/// residue before it. Proline has none, and neither has a chain's first
/// residue or one that follows a gap.
fn amide_hydrogen(backbones: &[Backbone], residue: usize) -> Option<[f32; 3]> {
    let backbone = backbones.get(residue)?;
    let previous = backbones.get(residue.checked_sub(1)?)?;
    if backbone.proline || !linked(previous, backbone) {
        return None;
    }
    let (carbon, oxygen, nitrogen) = (previous.carbon?, previous.oxygen?, backbone.nitrogen?);
    let axis = [
        carbon[0] - oxygen[0],
        carbon[1] - oxygen[1],
        carbon[2] - oxygen[2],
    ];
    let length = distance(carbon, oxygen);
    if length <= f32::EPSILON {
        return None;
    }
    let scale = NH_LENGTH / length;
    Some([
        nitrogen[0] + axis[0] * scale,
        nitrogen[1] + axis[1] * scale,
        nitrogen[2] + axis[2] * scale,
    ])
}

/// Whether `next` continues the chain from `previous` through a peptide bond.
fn linked(previous: &Backbone, next: &Backbone) -> bool {
    previous.chain == next.chain
        && match (previous.carbon, next.nitrogen) {
            (Some(carbon), Some(nitrogen)) => distance(carbon, nitrogen) <= PEPTIDE_BOND_LIMIT,
            _ => false,
        }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Ladder {
    Parallel,
    Antiparallel,
}

/// Per-residue flags for each pattern a residue can belong to.
struct Patterns {
    helix: Vec<bool>,
    three_ten: Vec<bool>,
    pi: Vec<bool>,
    turn: Vec<bool>,
    bridge: Vec<bool>,
    strand: Vec<bool>,
}

fn classify(
    backbones: &[Backbone],
    bonds: &Bonds,
    pairs: &[(usize, usize)],
    states: &mut [SecondaryStructure],
) {
    let mut patterns = Patterns {
        helix: vec![false; states.len()],
        three_ten: vec![false; states.len()],
        pi: vec![false; states.len()],
        turn: vec![false; states.len()],
        bridge: vec![false; states.len()],
        strand: vec![false; states.len()],
    };
    mark_turns_and_helices(bonds, &mut patterns);
    mark_bridges_and_ladders(bonds, pairs, &mut patterns);
    for (residue, state) in states.iter_mut().enumerate() {
        if *state == SecondaryStructure::Unknown {
            continue;
        }
        // DSSP's precedence: α-helix, isolated bridge, ladder strand, 3₁₀
        // helix, π-helix, turn, bend.
        *state = if patterns.helix[residue] {
            SecondaryStructure::AlphaHelix
        } else if patterns.bridge[residue] && !patterns.strand[residue] {
            SecondaryStructure::BetaBridge
        } else if patterns.strand[residue] {
            SecondaryStructure::Strand
        } else if patterns.three_ten[residue] {
            SecondaryStructure::ThreeTenHelix
        } else if patterns.pi[residue] {
            SecondaryStructure::PiHelix
        } else if patterns.turn[residue] {
            SecondaryStructure::Turn
        } else if bends(backbones, residue) {
            SecondaryStructure::Bend
        } else {
            *state
        };
    }
}

fn set(flags: &mut [bool], residues: impl IntoIterator<Item = usize>) {
    for residue in residues {
        if let Some(flag) = flags.get_mut(residue) {
            *flag = true;
        }
    }
}

/// Marks n-turns and the minimal helices two consecutive n-turns make.
fn mark_turns_and_helices(bonds: &Bonds, patterns: &mut Patterns) {
    let turn = |residue: usize, size: usize| bonds.contains(&(residue, residue + size));
    for residue in 0..patterns.turn.len() {
        for size in 3..=5 {
            if !turn(residue, size) {
                continue;
            }
            set(&mut patterns.turn, residue + 1..residue + size);
            // A minimal helix is two turns in a row; it covers the residues
            // between the first turn's carbonyl and the second turn's amide.
            if residue >= 1 && turn(residue - 1, size) {
                let target = match size {
                    3 => &mut patterns.three_ten,
                    4 => &mut patterns.helix,
                    _ => &mut patterns.pi,
                };
                set(target, residue..residue + size);
            }
        }
    }
}

/// Marks every bridged residue, and the residues of ladders: two bridges of
/// one kind in register.
fn mark_bridges_and_ladders(bonds: &Bonds, pairs: &[(usize, usize)], patterns: &mut Patterns) {
    let bond = |acceptor: Option<usize>, donor: Option<usize>| matches!((acceptor, donor), (Some(acceptor), Some(donor)) if bonds.contains(&(acceptor, donor)));
    let before = |residue: usize| residue.checked_sub(1);
    let mut bridges: IdentityHashMap<(usize, usize), Ladder> = IdentityHashMap::default();
    for &(i, j) in pairs {
        let parallel = (bond(before(i), Some(j)) && bond(Some(j), Some(i + 1)))
            || (bond(before(j), Some(i)) && bond(Some(i), Some(j + 1)));
        let antiparallel = (bond(Some(i), Some(j)) && bond(Some(j), Some(i)))
            || (bond(before(i), Some(j + 1)) && bond(before(j), Some(i + 1)));
        if parallel {
            bridges.insert((i, j), Ladder::Parallel);
        } else if antiparallel {
            bridges.insert((i, j), Ladder::Antiparallel);
        }
    }
    for (&(i, j), &ladder) in &bridges {
        set(&mut patterns.bridge, [i, j]);
        let next = match ladder {
            Ladder::Parallel => (i + 1, j + 1),
            Ladder::Antiparallel if j > 0 => (i + 1, j - 1),
            Ladder::Antiparallel => continue,
        };
        let key = (next.0.min(next.1), next.0.max(next.1));
        if bridges.get(&key) == Some(&ladder) {
            set(&mut patterns.strand, [i, j, next.0, next.1]);
        }
    }
}

/// The Cα-only fallback for a trace with no backbone atoms to bond: a helix
/// turns once in 3.6 residues, which puts Cα(i) and Cα(i+3) near 5 Å.
fn assign_zhang_skolnick(backbones: &[Backbone], states: &mut [SecondaryStructure]) {
    for index in 0..backbones.len() {
        let Some(window) = backbones.get(index..index + 4) else {
            break;
        };
        if window
            .iter()
            .any(|backbone| backbone.chain != window[0].chain)
        {
            continue;
        }
        if squared_distance(window[0].ca, window[3].ca) < 6.5 * 6.5
            && let Some(state) = states.get_mut(index..=index + 3)
        {
            for value in state {
                *value = SecondaryStructure::AlphaHelix;
            }
        }
    }
}

/// Whether the chain bends at `residue`: the Cα(i−2)→Cα(i) and Cα(i)→Cα(i+2)
/// directions differ by more than 70°, which is a Cα(i−2)–Cα(i)–Cα(i+2) angle
/// below 110°. All five residues must be linked without a break.
fn bends(backbones: &[Backbone], residue: usize) -> bool {
    let (Some(first), Some(last)) = (residue.checked_sub(2), residue.checked_add(2)) else {
        return false;
    };
    let Some(window) = backbones.get(first..=last) else {
        return false;
    };
    if !window.windows(2).all(|pair| linked(&pair[0], &pair[1])) {
        return false;
    }
    let (Some(before), Some(centre), Some(after)) = (window[0].ca, window[2].ca, window[4].ca)
    else {
        return false;
    };
    let incoming: [f64; 3] = core::array::from_fn(|axis| f64::from(centre[axis] - before[axis]));
    let outgoing: [f64; 3] = core::array::from_fn(|axis| f64::from(after[axis] - centre[axis]));
    let dot: f64 = (0..3).map(|axis| incoming[axis] * outgoing[axis]).sum();
    let norms = (0..3)
        .map(|axis| incoming[axis] * incoming[axis])
        .sum::<f64>()
        .sqrt()
        * (0..3)
            .map(|axis| outgoing[axis] * outgoing[axis])
            .sum::<f64>()
            .sqrt();
    norms > 0.0 && dot < norms * BEND_KAPPA_DEGREES.to_radians().cos()
}

fn squared_distance(left: Option<[f32; 3]>, right: Option<[f32; 3]>) -> f32 {
    let (Some(left), Some(right)) = (left, right) else {
        return f32::INFINITY;
    };
    left.iter().zip(right).map(|(a, b)| (a - b) * (a - b)).sum()
}

fn distance(left: [f32; 3], right: [f32; 3]) -> f32 {
    squared_distance(Some(left), Some(right)).sqrt()
}

#[cfg(test)]
#[path = "secondary_tests.rs"]
mod tests;
