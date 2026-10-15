//! The heavy-atom graph of a molecule and the ring structure read from it.
//!
//! Bridge finding is linear in bonds; ring extraction walks only the
//! non-bridge bonds, so acyclic molecules cost one pass.

use crate::aromaticity::{self, Aromaticity};
use crate::mol::MolRecord;
use molframe_core::BondOrder;

pub(super) struct Graph {
    pub(super) z: Vec<u8>,
    pub(super) charge: Vec<i8>,
    /// Heavy-atom bonds: (a, b, order 1..=4) with indices into the heavy atoms.
    pub(super) bonds: Vec<(usize, usize, u8)>,
    pub(super) adjacency: Vec<Vec<(usize, usize)>>,
    /// Hydrogens per heavy atom, when they can be counted.
    pub(super) hydrogens: Option<Vec<u8>>,
}

pub(super) fn graph(record: &MolRecord) -> Graph {
    let atoms = &record.molecule.atoms;
    let mut index = vec![usize::MAX; atoms.len()];
    let (mut z, mut charge) = (Vec::new(), Vec::new());
    for (i, atom) in atoms.iter().enumerate() {
        let n = atom.element.atomic_number();
        if n != 1 {
            index[i] = z.len();
            z.push(n);
            charge.push(
                match record.atom_metadata.get(i).and_then(|m| m.formal_charge) {
                    Some(formal) => formal,
                    None => 0,
                },
            );
        }
    }
    let mut explicit_h = vec![0_u8; z.len()];
    let mut any_h = false;
    let mut bonds = Vec::new();
    for b in &record.molecule.bonds {
        let (ia, ib) = (index[b.first], index[b.second]);
        match (ia != usize::MAX, ib != usize::MAX) {
            (true, true) => bonds.push((ia, ib, b.order)),
            (true, false) => {
                explicit_h[ia] = explicit_h[ia].saturating_add(1);
                any_h = true;
            }
            (false, true) => {
                explicit_h[ib] = explicit_h[ib].saturating_add(1);
                any_h = true;
            }
            (false, false) => {}
        }
    }
    let mut adjacency = vec![Vec::new(); z.len()];
    for (e, &(a, b, _)) in bonds.iter().enumerate() {
        adjacency[a].push((b, e));
        adjacency[b].push((a, e));
    }
    let aromatic_typed = bonds.iter().any(|b| b.2 == 4);
    let hydrogens = if any_h {
        Some(explicit_h)
    } else if aromatic_typed {
        None
    } else {
        implicit_hydrogens(&z, &charge, &bonds, z.len())
    };
    Graph {
        z,
        charge,
        bonds,
        adjacency,
        hydrogens,
    }
}

/// Hydrogens from standard valences; `None` when an element has no standard
/// valence here.
pub(super) fn implicit_hydrogens(
    z: &[u8],
    charge: &[i8],
    bonds: &[(usize, usize, u8)],
    n: usize,
) -> Option<Vec<u8>> {
    let mut sum = vec![0_i32; n];
    for &(a, b, o) in bonds {
        sum[a] += i32::from(o);
        sum[b] += i32::from(o);
    }
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let c = i32::from(charge[i]);
        let valences: &[i32] = match z[i] {
            6 => &[4],
            5 | 7 => &[3],
            8 => &[2],
            9 | 17 | 35 | 53 => &[1],
            16 => &[2, 4, 6],
            15 => &[3, 5],
            _ => return None,
        };
        // A charge moves the valence toward the neighbouring group: N+ and O+
        // take one more, N- and O- one fewer; C+ and C- one fewer.
        let shift = match z[i] {
            7 | 8 => c,
            6 => -c.abs(),
            _ => 0,
        };
        let valence = match valences.iter().map(|v| v + shift).find(|v| *v >= sum[i]) {
            Some(valence) => valence,
            None => sum[i],
        };
        out.push(match u8::try_from((valence - sum[i]).max(0)) {
            Ok(count) => count,
            Err(_) => 0,
        });
    }
    Some(out)
}

/// Perceives the aromatic atoms and bonds of the heavy-atom graph of `record`.
pub(super) fn perceive_aromaticity(record: &MolRecord, g: &Graph) -> Aromaticity {
    aromaticity::perceive(
        &record
            .molecule
            .atoms
            .iter()
            .filter(|a| a.element.atomic_number() != 1)
            .map(|a| a.element)
            .collect::<Vec<_>>(),
        &g.charge.iter().map(|c| Some(*c)).collect::<Vec<_>>(),
        &g.bonds
            .iter()
            .map(|&(a, b, o)| (a, b, bond_order(o)))
            .collect::<Vec<_>>(),
    )
}

fn bond_order(order: u8) -> BondOrder {
    match order {
        2 => BondOrder::Double,
        3 => BondOrder::Triple,
        4 => BondOrder::Aromatic,
        _ => BondOrder::Single,
    }
}

/// Bridges of the graph (bonds in no ring), by a depth-first low-link walk.
pub(super) fn bridges(g: &Graph) -> Vec<bool> {
    let n = g.z.len();
    let mut bridge = vec![false; g.bonds.len()];
    let (mut disc, mut low, mut timer) = (vec![0_usize; n], vec![0_usize; n], 1);
    for root in 0..n {
        if disc[root] != 0 {
            continue;
        }
        // (vertex, parent edge, next neighbour to visit)
        let mut stack = vec![(root, usize::MAX, 0_usize)];
        disc[root] = timer;
        low[root] = timer;
        timer += 1;
        while let Some(&mut (v, parent, ref mut next)) = stack.last_mut() {
            if *next < g.adjacency[v].len() {
                let (w, e) = g.adjacency[v][*next];
                *next += 1;
                if e == parent {
                    continue;
                }
                if disc[w] == 0 {
                    disc[w] = timer;
                    low[w] = timer;
                    timer += 1;
                    stack.push((w, e, 0));
                } else {
                    low[v] = low[v].min(disc[w]);
                }
            } else {
                stack.pop();
                if let Some(&(p, _, _)) = stack.last() {
                    low[p] = low[p].min(low[v]);
                    if low[v] > disc[p] {
                        bridge[parent] = true;
                    }
                }
            }
        }
    }
    bridge
}

/// A smallest set of smallest rings, each as the bonds it uses.
pub(super) fn rings(g: &Graph, bridge: &[bool]) -> Vec<Vec<usize>> {
    let n = g.z.len();
    // Cycle rank: bonds - atoms + connected components.
    let mut seen = vec![false; n];
    let mut components = 0;
    for s in 0..n {
        if seen[s] {
            continue;
        }
        components += 1;
        let mut stack = vec![s];
        seen[s] = true;
        while let Some(v) = stack.pop() {
            for &(w, _) in &g.adjacency[v] {
                if !seen[w] {
                    seen[w] = true;
                    stack.push(w);
                }
            }
        }
    }
    let rank = (g.bonds.len() + components).saturating_sub(n);
    if rank == 0 {
        return Vec::new();
    }
    // Candidates: for every ring bond, the shortest cycle through it.
    let mut candidates: Vec<Vec<usize>> = Vec::new();
    for (e, &(a, b, _)) in g.bonds.iter().enumerate() {
        if bridge[e] {
            continue;
        }
        let mut prev = vec![(usize::MAX, usize::MAX); n];
        let mut queue = std::collections::VecDeque::from([a]);
        prev[a] = (a, usize::MAX);
        while let Some(v) = queue.pop_front() {
            if v == b {
                break;
            }
            for &(w, f) in &g.adjacency[v] {
                if f == e || prev[w].0 != usize::MAX {
                    continue;
                }
                prev[w] = (v, f);
                queue.push_back(w);
            }
        }
        if prev[b].0 == usize::MAX {
            continue;
        }
        let mut cycle = vec![e];
        let mut at = b;
        while at != a {
            cycle.push(prev[at].1);
            at = prev[at].0;
        }
        cycle.sort_unstable();
        if !candidates.contains(&cycle) {
            candidates.push(cycle);
        }
    }
    candidates.sort_by_key(Vec::len);
    // Keep the smallest that are independent over GF(2).
    let words = g.bonds.len().div_ceil(64);
    let mut basis: Vec<Vec<u64>> = Vec::new();
    let mut chosen = Vec::new();
    for cycle in candidates {
        let mut v = vec![0_u64; words];
        for &b in &cycle {
            v[b / 64] ^= 1 << (b % 64);
        }
        for row in &basis {
            let lead = row
                .iter()
                .position(|w| *w != 0)
                .map(|i| (i, row[i].trailing_zeros()));
            if let Some((i, bit)) = lead
                && v[i] >> bit & 1 == 1
            {
                for (x, y) in v.iter_mut().zip(row) {
                    *x ^= y;
                }
            }
        }
        if v.iter().any(|w| *w != 0) {
            basis.push(v);
            chosen.push(cycle);
            if chosen.len() == rank {
                break;
            }
        }
    }
    chosen
}

/// A bond that is in a ring is not a bridge.
pub(super) const fn bridge_is_ring(bridge: &[bool], e: usize) -> bool {
    !bridge[e]
}
