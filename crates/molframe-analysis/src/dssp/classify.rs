//! Naming residues from the backbone hydrogen-bond pattern.
//!
//! The rules follow Kabsch and Sander. A single `i→i+helix_offset` bond
//! brackets an α-helix, as this kernel has always read it. A 3₁₀ or π helix
//! needs two consecutive turns of its size, because a lone short bond is a turn
//! rather than a helix. Reciprocal or offset bonds between distant residues are
//! β-bridges; two bridges in register are a ladder, whose residues are strand,
//! and a bridge with no partner in register is an isolated β-bridge. What
//! remains is a turn when a short bond brackets it, a bend when the Cα trace
//! turns sharply there, and coil otherwise. Where classes overlap the DSSP
//! precedence decides: α-helix, β-bridge, strand, 3₁₀, π, turn, bend.

use super::DsspOptions;
use molframe_core::SecondaryStructure;
use std::collections::BTreeSet;

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Ladder {
    Parallel,
    Antiparallel,
}

/// Assigns a state to each residue from the hydrogen-bond pattern.
///
/// `alpha_carbons` carries each residue's Cα, used only for bends; a residue
/// with no Cα, or with a neighbour two away missing one, never bends.
pub(super) fn classify(
    bonds: &BTreeSet<(usize, usize)>,
    evaluable: &[bool],
    alpha_carbons: &[Option<[f32; 3]>],
    options: &DsspOptions,
) -> Vec<SecondaryStructure> {
    let count = evaluable.len();
    let has = |i: usize, j: usize| bonds.contains(&(i, j));
    let mut alpha = vec![false; count];
    let mut three_ten = vec![false; count];
    let mut pi = vec![false; count];
    let mut turn = vec![false; count];

    for &(i, j) in bonds {
        if j <= i {
            continue;
        }
        let size = j - i;
        if size == options.helix_offset {
            fill(&mut alpha, i + 1..j);
        }
        if options.turn_offsets.contains(&size) {
            fill(&mut turn, i + 1..j);
        }
        // Two consecutive n-turns at i−1 and i make a minimal n-helix over
        // residues i..i+n−1.
        if i >= 1 && has(i - 1, j - 1) {
            if size == options.three_ten_offset {
                fill(&mut three_ten, i..j);
            }
            if size == options.pi_offset {
                fill(&mut pi, i..j);
            }
        }
    }

    let bridges = bridges(bonds, count);
    let mut paired = vec![false; count];
    let mut ladder = vec![false; count];
    for &(i, j, kind) in &bridges {
        paired[i] = true;
        paired[j] = true;
        let in_register = match kind {
            Ladder::Parallel => [
                (i.checked_add(1), j.checked_add(1)),
                (i.checked_sub(1), j.checked_sub(1)),
            ],
            Ladder::Antiparallel => [
                (i.checked_add(1), j.checked_sub(1)),
                (i.checked_sub(1), j.checked_add(1)),
            ],
        };
        let partnered = in_register.iter().any(|pair| match *pair {
            (Some(a), Some(b)) => bridges.contains(&(a.min(b), a.max(b), kind)),
            _ => false,
        });
        if partnered {
            ladder[i] = true;
            ladder[j] = true;
        }
    }

    (0..count)
        .map(|residue| {
            if !evaluable[residue] {
                SecondaryStructure::Unknown
            } else if alpha[residue] {
                SecondaryStructure::AlphaHelix
            } else if paired[residue] && !ladder[residue] {
                SecondaryStructure::BetaBridge
            } else if ladder[residue] {
                SecondaryStructure::Strand
            } else if three_ten[residue] {
                SecondaryStructure::ThreeTenHelix
            } else if pi[residue] {
                SecondaryStructure::PiHelix
            } else if turn[residue] {
                SecondaryStructure::Turn
            } else if bends(alpha_carbons, residue, options.bend_angle_degrees) {
                SecondaryStructure::Bend
            } else {
                SecondaryStructure::Coil
            }
        })
        .collect()
}

fn fill(flags: &mut [bool], range: std::ops::Range<usize>) {
    let end = range.end.min(flags.len());
    for flag in flags.iter_mut().take(end).skip(range.start) {
        *flag = true;
    }
}

/// Every β-bridge as `(lower, higher, kind)`, from the bonds between residues
/// more than two apart.
fn bridges(bonds: &BTreeSet<(usize, usize)>, count: usize) -> BTreeSet<(usize, usize, Ladder)> {
    let has = |i: usize, j: usize| bonds.contains(&(i, j));
    let mut found = BTreeSet::new();
    for &(a, b) in bonds {
        let (i, j) = (a.min(b), a.max(b));
        if j - i <= 2 {
            continue;
        }
        let antiparallel = (has(i, j) && has(j, i))
            || (i >= 1 && j + 1 < count && has(i - 1, j + 1) && has(j - 1, i + 1));
        let parallel = (i >= 1 && has(i - 1, j) && has(j, i + 1))
            || (j >= 1 && has(j - 1, i) && has(i, j + 1));
        if parallel {
            found.insert((i, j, Ladder::Parallel));
        } else if antiparallel {
            found.insert((i, j, Ladder::Antiparallel));
        }
    }
    found
}

/// Whether the Cα direction changes by more than `threshold` degrees at
/// `residue`, measured between Cα(i−2)→Cα(i) and Cα(i)→Cα(i+2).
fn bends(alpha_carbons: &[Option<[f32; 3]>], residue: usize, threshold: f32) -> bool {
    let (Some(first), Some(last)) = (residue.checked_sub(2), residue.checked_add(2)) else {
        return false;
    };
    let (Some(Some(before)), Some(Some(centre)), Some(Some(after))) = (
        alpha_carbons.get(first),
        alpha_carbons.get(residue),
        alpha_carbons.get(last),
    ) else {
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
    norms > 0.0 && dot < norms * f64::from(threshold).to_radians().cos()
}
