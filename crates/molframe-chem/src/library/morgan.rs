//! Morgan-style fingerprints, Tanimoto similarity and Butina clustering.
//!
//! Fingerprint construction is linear in atoms times radius; clustering builds
//! the full neighbour table, so it is quadratic in the number of molecules.

use super::graph::{bridges, graph, perceive_aromaticity};
use crate::mol::MolRecord;

/// Bits in a [`Fingerprint`].
pub const FINGERPRINT_BITS: usize = 2048;
const BITS: u64 = 2048;
const WORDS: usize = FINGERPRINT_BITS / 64;

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
    let arom = perceive_aromaticity(record, &g);
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
        let bit = match usize::try_from(id % BITS) {
            Ok(bit) => bit,
            Err(_) => 0,
        };
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
        let densest = left.iter().max_by_key(|&&i| {
            (
                neighbours[i].iter().filter(|&&j| !taken[j]).count(),
                std::cmp::Reverse(i),
            )
        });
        let centre = match densest {
            Some(centre) => *centre,
            None => left[0],
        };
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
