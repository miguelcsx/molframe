//! Chemistry-derived side-chain torsion paths.

use crate::Component;
use pdbiox_core::{BondOrder, Element};

const MAX_PATH_ATOMS: usize = 8;

/// Ordered component atoms defining up to χ1–χ5.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SideChainDefinition {
    /// `N`, `CA`, then the selected heavy-atom path from `CB` outward.
    pub atoms: Box<[Box<str>]>,
}

impl SideChainDefinition {
    /// Number of consecutive χ torsions defined by this path.
    #[must_use]
    pub fn torsion_count(&self) -> usize {
        self.atoms.len().saturating_sub(3)
    }
}

/// Derives the conventional deterministic χ path from CCD connectivity.
///
/// Branch ties use component atom identifiers, so equivalent terminal atoms
/// never depend on row order. Delocalised/multiple bonds may terminate a path
/// but remain its last atom, which retains χ2 for aromatic residues and χ5 for
/// arginine-like guanidinium groups.
#[must_use]
pub fn side_chain_definition(component: &Component) -> Option<SideChainDefinition> {
    let nitrogen = atom_index(component, "N")?;
    let alpha = atom_index(component, "CA")?;
    let beta = atom_index(component, "CB")?;
    let adjacency = adjacency(component);
    if !connected(&adjacency, nitrogen, alpha) || !connected(&adjacency, alpha, beta) {
        return None;
    }
    let mut path = vec![nitrogen, alpha, beta];
    let mut best = path.clone();
    extend(component, &adjacency, &mut path, &mut best);
    Some(SideChainDefinition {
        atoms: best
            .into_iter()
            .filter_map(|index| component.atoms.get(index).map(|atom| atom.name.clone()))
            .collect(),
    })
}

fn extend(
    component: &Component,
    adjacency: &[Vec<(usize, BondOrder)>],
    path: &mut Vec<usize>,
    best: &mut Vec<usize>,
) {
    if better(component, path, best) {
        best.clone_from(path);
    }
    if path.len() >= MAX_PATH_ATOMS {
        return;
    }
    let Some((&previous, prefix)) = path.split_last() else {
        return;
    };
    let Some(&before) = prefix.last() else {
        return;
    };
    let incoming = adjacency
        .get(before)
        .and_then(|neighbors| neighbors.iter().find(|(atom, _)| *atom == previous))
        .map(|(_, order)| *order);
    if !incoming.is_some_and(rotatable) {
        return;
    }
    let mut candidates = match adjacency.get(previous) {
        Some(candidates) => candidates.clone(),
        None => Vec::new(),
    };
    candidates.sort_by(|(left, _), (right, _)| {
        atom_name(component, *left).cmp(atom_name(component, *right))
    });
    for (next, _) in candidates {
        if path.contains(&next) || !side_chain_atom(component, next) {
            continue;
        }
        path.push(next);
        extend(component, adjacency, path, best);
        path.pop();
    }
}

fn adjacency(component: &Component) -> Vec<Vec<(usize, BondOrder)>> {
    let mut output = vec![Vec::new(); component.atoms.len()];
    for bond in component.bonds.iter() {
        let (Some(left), Some(right)) = (
            atom_index(component, &bond.atom_a),
            atom_index(component, &bond.atom_b),
        ) else {
            continue;
        };
        output[left].push((right, bond.order));
        output[right].push((left, bond.order));
    }
    output
}

fn atom_index(component: &Component, name: &str) -> Option<usize> {
    component
        .atoms
        .iter()
        .position(|atom| atom.name.as_ref() == name)
}

fn atom_name(component: &Component, index: usize) -> &str {
    component
        .atoms
        .get(index)
        .map_or("", |atom| atom.name.as_ref())
}

fn side_chain_atom(component: &Component, index: usize) -> bool {
    component.atoms.get(index).is_some_and(|atom| {
        atom.element != Element::HYDROGEN
            && !matches!(atom.name.as_ref(), "N" | "CA" | "C" | "O" | "OXT")
    })
}

fn connected(adjacency: &[Vec<(usize, BondOrder)>], left: usize, right: usize) -> bool {
    adjacency
        .get(left)
        .is_some_and(|neighbors| neighbors.iter().any(|(atom, _)| *atom == right))
}

fn rotatable(order: BondOrder) -> bool {
    matches!(order, BondOrder::Single | BondOrder::Unknown)
}

fn better(component: &Component, candidate: &[usize], current: &[usize]) -> bool {
    candidate.len() > current.len()
        || (candidate.len() == current.len()
            && candidate
                .iter()
                .map(|index| atom_name(component, *index))
                .cmp(current.iter().map(|index| atom_name(component, *index)))
                .is_lt())
}

#[cfg(test)]
#[path = "side_chain_tests.rs"]
mod tests;
