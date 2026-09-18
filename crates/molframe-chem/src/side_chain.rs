//! Chemistry-derived side-chain torsion paths.

use crate::Component;
use molframe_core::{BondOrder, Element};

const MAX_PATH_ATOMS: usize = 8;

/// Ordered component atoms defining up to χ1–χ5.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SideChainDefinition {
    /// Annotated nitrogen and alpha-carbon anchors, then the side-chain path.
    pub atoms: Box<[Box<str>]>,
}

/// Component-local atom roles resolved by a CCD-aware or caller-supplied profile.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SideChainRoles<'a> {
    /// Protein peptide nitrogen atom.
    pub nitrogen: &'a str,
    /// Protein alpha-carbon atom.
    pub alpha_carbon: &'a str,
    /// Atoms belonging to the side chain, including its first heavy atom.
    pub side_chain_atoms: &'a [&'a str],
}

impl SideChainDefinition {
    /// Number of consecutive χ torsions defined by this path.
    #[must_use]
    pub fn torsion_count(&self) -> usize {
        if self.atoms.len() < 3 {
            0
        } else {
            self.atoms.len() - 3
        }
    }
}

/// Derives a deterministic χ path from CCD connectivity and explicit roles.
///
/// Branch ties use component atom identifiers, so equivalent terminal atoms
/// never depend on row order. Delocalised/multiple bonds may terminate a path
/// but remain its last atom, which retains χ2 for aromatic residues and χ5 for
/// arginine-like guanidinium groups.
#[must_use]
pub fn side_chain_definition(
    component: &Component,
    roles: &SideChainRoles<'_>,
) -> Option<SideChainDefinition> {
    let nitrogen = atom_index(component, roles.nitrogen)?;
    let alpha = atom_index(component, roles.alpha_carbon)?;
    let adjacency = adjacency(component);
    if !connected(&adjacency, nitrogen, alpha) {
        return None;
    }
    let beta = first_side_chain_neighbour(component, &adjacency, alpha, roles)?;
    let mut path = vec![nitrogen, alpha, beta];
    let mut best = path.clone();
    extend(component, roles, &adjacency, &mut path, &mut best);
    Some(SideChainDefinition {
        atoms: best
            .into_iter()
            .filter_map(|index| component.atoms.get(index).map(|atom| atom.name.clone()))
            .collect(),
    })
}

fn extend(
    component: &Component,
    roles: &SideChainRoles<'_>,
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
        if path.contains(&next) || !side_chain_atom(component, roles, next) {
            continue;
        }
        path.push(next);
        extend(component, roles, adjacency, path, best);
        path.pop();
    }
}

fn first_side_chain_neighbour(
    component: &Component,
    adjacency: &[Vec<(usize, BondOrder)>],
    alpha: usize,
    roles: &SideChainRoles<'_>,
) -> Option<usize> {
    adjacency
        .get(alpha)?
        .iter()
        .map(|(index, _)| *index)
        .filter(|index| side_chain_atom(component, roles, *index))
        .min_by(|left, right| atom_name(component, *left).cmp(atom_name(component, *right)))
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

fn side_chain_atom(component: &Component, roles: &SideChainRoles<'_>, index: usize) -> bool {
    component.atoms.get(index).is_some_and(|atom| {
        atom.element != Element::HYDROGEN && roles.side_chain_atoms.contains(&atom.name.as_ref())
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
