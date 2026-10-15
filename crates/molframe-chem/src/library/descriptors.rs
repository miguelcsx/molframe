//! Screening descriptors for a small molecule and the rule-of-five tally.
//!
//! Everything is a graph computation over a [`MolRecord`]; each result states
//! what it rests on, and a quantity that cannot be established is absent
//! rather than guessed.

use super::graph::{bridge_is_ring, bridges, graph, perceive_aromaticity, rings};
use crate::mol::MolRecord;

/// Mass of a hydrogen, standard atomic weight.
const HYDROGEN_WEIGHT: f64 = 1.008;

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

/// Describes `record`.
#[must_use]
pub fn describe(record: &MolRecord) -> LigandDescriptors {
    let g = graph(record);
    let n = g.z.len();
    let bridge = bridges(&g);
    let ring_set = rings(&g, &bridge);
    let arom = perceive_aromaticity(record, &g);
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
    let small = |count: usize| match u32::try_from(count) {
        Ok(count) => count,
        Err(_) => u32::MAX,
    };
    f64::from(small(a)) / f64::from(small(b))
}
