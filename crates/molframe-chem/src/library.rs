//! Small-molecule library chemistry: the properties a screening campaign filters
//! on, a Morgan-style fingerprint, and similarity clustering.
//!
//! Everything here is a graph computation over a [`MolRecord`]; no geometry is
//! consulted. Each result says what it did not establish: hydrogens are taken
//! from the file when it has explicit ones and otherwise from standard valences
//! (which is not possible for aromatic-typed bonds without explicit hydrogens,
//! so those molecules get no mass or donor count rather than a guess), the
//! aromatic ring count is the conservative one of `molframe-chem`'s own
//! aromaticity perception, and no octanol-water partition is estimated.

use crate::aromaticity;
use crate::mol::MolRecord;
use molframe_core::BondOrder;

/// Mass of a hydrogen, standard atomic weight.
const HYDROGEN_WEIGHT: f64 = 1.008;
/// Bits in a [`Fingerprint`].
pub const FINGERPRINT_BITS: usize = 2048;
const BITS: u64 = 2048;
const WORDS: usize = FINGERPRINT_BITS / 64;

/// What a screening campaign filters on, with what each number rests on.
#[derive(Clone, Debug, PartialEq)]
pub struct LigandDescriptors {
    /// Atoms other than hydrogen.
    pub heavy_atoms: usize,
    /// Molecular weight, g/mol; absent when the hydrogens cannot be counted.
    pub molecular_weight: Option<f64>,
    /// N and O atoms carrying a hydrogen (the rule of five's own count).
    pub hydrogen_bond_donors: Option<usize>,
    /// N and O atoms (the rule of five's own count).
    pub hydrogen_bond_acceptors: usize,
    /// Single, acyclic bonds between two atoms that each have another heavy
    /// neighbour, not touching a triple bond.
    pub rotatable_bonds: usize,
    /// Rings in a smallest set of smallest rings.
    pub rings: usize,
    /// Rings whose every bond is aromatic by the conservative graph rule: a
    /// lower bound (ring systems with exocyclic multiple bonds are not inferred).
    pub aromatic_rings: usize,
    /// Fraction of carbons that are sp3 (single bonds only, not aromatic).
    pub fraction_csp3: Option<f64>,
    /// Sum of the declared formal charges.
    pub formal_charge: i32,
}

/// Lipinski's rule of five, as far as it can be checked here.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RuleOfFive {
    /// Rules broken among those checked.
    pub violations: usize,
    /// Rules that could be checked: mass, donors and acceptors. The fourth
    /// (partition coefficient) is not estimated, so it is never counted.
    pub checked: usize,
}

impl LigandDescriptors {
    /// The rule of five over what is known.
    #[must_use]
    pub fn rule_of_five(&self) -> RuleOfFive {
        let mut violations = 0;
        let mut checked = 1; // acceptors always
        if self.hydrogen_bond_acceptors > 10 {
            violations += 1;
        }
        if let Some(mw) = self.molecular_weight {
            checked += 1;
            if mw > 500.0 {
                violations += 1;
            }
        }
        if let Some(d) = self.hydrogen_bond_donors {
            checked += 1;
            if d > 5 {
                violations += 1;
            }
        }
        RuleOfFive {
            violations,
            checked,
        }
    }
}

struct Graph {
    z: Vec<u8>,
    charge: Vec<i8>,
    /// Heavy-atom bonds: (a, b, order 1..=4) with indices into the heavy atoms.
    bonds: Vec<(usize, usize, u8)>,
    adjacency: Vec<Vec<(usize, usize)>>,
    /// Hydrogens per heavy atom, when they can be counted.
    hydrogens: Option<Vec<u8>>,
}

fn graph(record: &MolRecord) -> Graph {
    let atoms = &record.molecule.atoms;
    let mut index = vec![usize::MAX; atoms.len()];
    let (mut z, mut charge) = (Vec::new(), Vec::new());
    for (i, atom) in atoms.iter().enumerate() {
        let n = atom.element.atomic_number();
        if n != 1 {
            index[i] = z.len();
            z.push(n);
            charge.push(
                record
                    .atom_metadata
                    .get(i)
                    .and_then(|m| m.formal_charge)
                    .unwrap_or(0),
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
fn implicit_hydrogens(
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
        let valence = valences
            .iter()
            .map(|v| v + shift)
            .find(|v| *v >= sum[i])
            .unwrap_or(sum[i]);
        out.push(u8::try_from((valence - sum[i]).max(0)).unwrap_or(0));
    }
    Some(out)
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
fn bridges(g: &Graph) -> Vec<bool> {
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
fn rings(g: &Graph, bridge: &[bool]) -> Vec<Vec<usize>> {
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

/// Describes `record`.
#[must_use]
pub fn describe(record: &MolRecord) -> LigandDescriptors {
    let g = graph(record);
    let n = g.z.len();
    let bridge = bridges(&g);
    let ring_set = rings(&g, &bridge);
    let arom = aromaticity::perceive(
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
    );
    let aromatic_rings = ring_set
        .iter()
        .filter(|r| r.iter().all(|b| arom.bonds[*b]))
        .count();
    let weights = |z: u8| {
        crate::element_properties(molframe_core::Element::from_atomic_number(z))
            .map(|p| p.atomic_weight)
    };
    let molecular_weight = g.hydrogens.as_ref().and_then(|h| {
        let mut total = 0.0;
        for (z, hydrogens) in g.z.iter().zip(h) {
            total += weights(*z)? + f64::from(*hydrogens) * HYDROGEN_WEIGHT;
        }
        Some(total)
    });
    let donors = g.hydrogens.as_ref().map(|h| {
        (0..n)
            .filter(|i| matches!(g.z[*i], 7 | 8) && h[*i] > 0)
            .count()
    });
    let acceptors = g.z.iter().filter(|z| matches!(**z, 7 | 8)).count();
    let triple: Vec<bool> = (0..n)
        .map(|i| g.adjacency[i].iter().any(|&(_, e)| g.bonds[e].2 == 3))
        .collect();
    let rotatable = g
        .bonds
        .iter()
        .enumerate()
        .filter(|&(e, &(a, b, o))| {
            o == 1
                && !bridge_is_ring(&bridge, e)
                && g.adjacency[a].len() > 1
                && g.adjacency[b].len() > 1
                && !triple[a]
                && !triple[b]
        })
        .count();
    let carbons: Vec<usize> = (0..n).filter(|i| g.z[*i] == 6).collect();
    let fsp3 = (!carbons.is_empty()).then(|| {
        let sp3 = carbons
            .iter()
            .filter(|&&i| !arom.atoms[i] && g.adjacency[i].iter().all(|&(_, e)| g.bonds[e].2 == 1))
            .count();
        ratio(sp3, carbons.len())
    });
    LigandDescriptors {
        heavy_atoms: n,
        molecular_weight,
        hydrogen_bond_donors: donors,
        hydrogen_bond_acceptors: acceptors,
        rotatable_bonds: rotatable,
        rings: ring_set.len(),
        aromatic_rings,
        fraction_csp3: fsp3,
        formal_charge: g.charge.iter().map(|c| i32::from(*c)).sum(),
    }
}

/// `a / b` for counts far below 2^32.
fn ratio(a: usize, b: usize) -> f64 {
    f64::from(u32::try_from(a).unwrap_or(u32::MAX))
        / f64::from(u32::try_from(b).unwrap_or(u32::MAX))
}

/// A bond that is in a ring is not a bridge.
const fn bridge_is_ring(bridge: &[bool], e: usize) -> bool {
    !bridge[e]
}

/// A fixed-length bit vector of the environments in a molecule.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Fingerprint([u64; WORDS]);

impl Fingerprint {
    /// Bits set.
    #[must_use]
    pub fn count(&self) -> u32 {
        self.0.iter().map(|w| w.count_ones()).sum()
    }
}

fn mix(h: u64, v: u64) -> u64 {
    // FNV-1a over the eight bytes of v.
    let mut h = h;
    for byte in v.to_le_bytes() {
        h ^= u64::from(byte);
        h = h.wrapping_mul(0x0000_0100_0000_01B3);
    }
    h
}

/// A Morgan-style (circular) fingerprint of radius 2, hashed to
/// [`FINGERPRINT_BITS`] bits: each atom's environment out to two bonds sets one
/// bit. It follows the construction of extended-connectivity fingerprints but
/// does not reproduce any toolkit's bits, so compare fingerprints made here with
/// each other only.
#[must_use]
pub fn fingerprint(record: &MolRecord) -> Fingerprint {
    let g = graph(record);
    let n = g.z.len();
    let bridge = bridges(&g);
    let arom = aromaticity::perceive(
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
    );
    let in_ring: Vec<bool> = (0..n)
        .map(|i| g.adjacency[i].iter().any(|&(_, e)| !bridge[e]))
        .collect();
    let mut ids: Vec<u64> = (0..n)
        .map(|i| {
            let h = g.hydrogens.as_ref().map_or(255, |h| u64::from(h[i]));
            let mut x = mix(0xcbf2_9ce4_8422_2325, u64::from(g.z[i]));
            x = mix(x, g.adjacency[i].len() as u64);
            x = mix(x, h);
            x = mix(x, u64::from(g.charge[i].unsigned_abs()));
            x = mix(x, u64::from(g.charge[i] < 0));
            x = mix(x, u64::from(in_ring[i]));
            mix(x, u64::from(arom.atoms[i]))
        })
        .collect();
    let mut bits = [0_u64; WORDS];
    let mut set = |id: u64| {
        let bit = usize::try_from(id % BITS).unwrap_or(0);
        bits[bit / 64] |= 1 << (bit % 64);
    };
    for &id in &ids {
        set(id);
    }
    for _ in 0..2 {
        let next: Vec<u64> = (0..n)
            .map(|i| {
                let mut around: Vec<(u64, u64)> = g.adjacency[i]
                    .iter()
                    .map(|&(w, e)| {
                        let order = if arom.bonds[e] {
                            4
                        } else {
                            u64::from(g.bonds[e].2)
                        };
                        (order, ids[w])
                    })
                    .collect();
                around.sort_unstable();
                let mut x = mix(0xcbf2_9ce4_8422_2325, ids[i]);
                for (o, id) in around {
                    x = mix(mix(x, o), id);
                }
                x
            })
            .collect();
        ids = next;
        for &id in &ids {
            set(id);
        }
    }
    Fingerprint(bits)
}

/// Tanimoto similarity of two fingerprints: bits in common over bits in either.
/// Two empty fingerprints are not similar (0), not identical.
#[must_use]
pub fn tanimoto(a: &Fingerprint, b: &Fingerprint) -> f64 {
    let (mut both, mut either) = (0_u32, 0_u32);
    for (x, y) in a.0.iter().zip(&b.0) {
        both += (x & y).count_ones();
        either += (x | y).count_ones();
    }
    if either == 0 {
        0.0
    } else {
        f64::from(both) / f64::from(either)
    }
}

/// Butina clusters: repeatedly take the molecule with the most neighbours at
/// Tanimoto similarity ≥ `threshold` as a centre, with every unclustered
/// neighbour. Each cluster lists its centre first; every index is in exactly one.
#[must_use]
pub fn cluster(prints: &[Fingerprint], threshold: f64) -> Vec<Vec<usize>> {
    let n = prints.len();
    let neighbours: Vec<Vec<usize>> = (0..n)
        .map(|i| {
            (0..n)
                .filter(|&j| j != i && tanimoto(&prints[i], &prints[j]) >= threshold)
                .collect()
        })
        .collect();
    let mut left: Vec<usize> = (0..n).collect();
    let mut taken = vec![false; n];
    let mut out = Vec::new();
    while !left.is_empty() {
        let centre = *left
            .iter()
            .max_by_key(|&&i| {
                (
                    neighbours[i].iter().filter(|&&j| !taken[j]).count(),
                    std::cmp::Reverse(i),
                )
            })
            .unwrap_or(&left[0]);
        let mut members = vec![centre];
        taken[centre] = true;
        for &j in &neighbours[centre] {
            if !taken[j] {
                taken[j] = true;
                members.push(j);
            }
        }
        left.retain(|i| !taken[*i]);
        out.push(members);
    }
    out
}

#[cfg(test)]
#[path = "library_tests.rs"]
mod tests;
