//! Native DSSP 4 pattern assignment, including stretch-level helix precedence.

use super::geometry::{bends, continuous, pp_conformation};
use super::{Bonds, DsspBackbone as Backbone, DsspOptions, bridges};
use molframe_core::SecondaryStructure as Ss;

pub(super) fn classify(
    backbones: &[Backbone],
    bonds: &Bonds,
    pairs: &[(usize, usize)],
    states: &mut [Ss],
    options: &DsspOptions,
) {
    bridges::assign(backbones, bonds, pairs, states);
    let turn = |i: usize, size: usize| {
        i.checked_add(size)
            .is_some_and(|j| bonds.contains(&(i, j)) && continuous(backbones, i, j))
    };
    for (size, kind) in [
        (options.helix_offset, Ss::AlphaHelix),
        (options.three_ten_offset, Ss::ThreeTenHelix),
        (options.pi_offset, Ss::PiHelix),
    ] {
        for i in 1..states.len().saturating_sub(size) {
            if !turn(i - 1, size) || !turn(i, size) {
                continue;
            }
            let stretch = &mut states[i..i + size];
            let accepted = match kind {
                Ss::AlphaHelix => true,
                Ss::ThreeTenHelix => stretch
                    .iter()
                    .all(|s| matches!(s, Ss::Coil | Ss::ThreeTenHelix)),
                Ss::PiHelix => stretch
                    .iter()
                    .all(|s| matches!(s, Ss::Coil | Ss::PiHelix | Ss::AlphaHelix)),
                _ => false,
            };
            if accepted {
                stretch.fill(kind);
            }
        }
    }
    for i in 1..states.len().saturating_sub(1) {
        if states[i] != Ss::Coil {
            continue;
        }
        let is_turn = options
            .turn_offsets
            .clone()
            .any(|size| (1..size.min(i + 1)).any(|offset| turn(i - offset, size)));
        if is_turn {
            states[i] = Ss::Turn;
        } else if bends(backbones, i, options.bend_angle_degrees) {
            states[i] = Ss::Bend;
        }
    }
    // PPII is a phi/psi conformation, not a proline sequence test. It is the
    // final assignment pass and must not erase a hydrogen-bonded turn or bend.
    let qualified: Vec<bool> = (0..states.len())
        .map(|i| pp_conformation(backbones, i))
        .collect();
    for i in 1..states.len().saturating_sub(3) {
        if qualified[i..i + 3].iter().all(|q| *q) {
            for state in &mut states[i..i + 3] {
                if *state == Ss::Coil {
                    *state = Ss::PolyProline;
                }
            }
        }
    }
}
