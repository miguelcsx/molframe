//! Interchangeable atoms of a residue, from the bond graph of its component.
//!
//! Two atoms are interchangeable when some automorphism of the residue's
//! heavy-atom graph exchanges them. The graph is taken with bond orders
//! ignored: a dictionary draws a carboxylate with one double and one single
//! bond, so an order-aware graph would call its two oxygens different although
//! no experiment can tell them apart. Ignoring orders also exchanges the two
//! methyls of valine and leucine, which are distinguishable, so an automorphism
//! is kept only when every atom it moves hangs from a fixed planar centre, one
//! with a double or aromatic bond. Results are cached by component and by the
//! set of atoms present, because most residues of a structure repeat both.

use crate::{CompareError, component_fragment};
use molframe_chem::{Component, automorphisms};
use molframe_core::BondOrder;
use std::collections::BTreeMap;
use std::sync::Arc;

/// Non-identity exchanges over the present atoms of a component, indexed by
/// position in that component's present-atom order.
pub(super) type Permutations = Arc<[Box<[u32]>]>;

type CacheKey = (Box<str>, Box<[u32]>);

#[derive(Default)]
pub(super) struct SymmetryCache {
    entries: BTreeMap<CacheKey, Option<Permutations>>,
}

impl SymmetryCache {
    /// The admitted exchanges of the atoms at `present` positions of
    /// `component`, or `None` when enumeration exceeded `limit`.
    pub(super) fn permutations(
        &mut self,
        component: &Component,
        present: &[u32],
        limit: usize,
    ) -> Option<Permutations> {
        let key = (component.id.clone(), Box::from(present));
        self.entries
            .entry(key)
            .or_insert_with(|| admitted(component, present, limit).ok())
            .clone()
    }
}

fn admitted(
    component: &Component,
    present: &[u32],
    limit: usize,
) -> Result<Permutations, CompareError> {
    let wanted = |name: &str| {
        component
            .atoms
            .iter()
            .position(|atom| atom.name.as_ref() == name)
            .is_some_and(|position| present.iter().any(|kept| *kept as usize == position))
    };
    let fragment = component_fragment(component, wanted);
    let agnostic = Component {
        atoms: fragment
            .atoms
            .iter()
            .cloned()
            .map(|mut atom| {
                atom.aromatic = false;
                atom
            })
            .collect(),
        bonds: fragment
            .bonds
            .iter()
            .cloned()
            .map(|mut bond| {
                bond.order = BondOrder::Single;
                bond.aromatic = false;
                bond
            })
            .collect(),
        ..fragment.clone()
    };
    let all = automorphisms(&agnostic, limit)
        .map_err(|error| CompareError::MappingLimit { limit: error.limit })?;
    let guard = Guard::new(component, &fragment);
    let kept: Vec<Box<[u32]>> = all
        .into_iter()
        .filter(|mapping| !is_identity(mapping) && guard.admits(mapping))
        .collect();
    Ok(kept.into())
}

fn is_identity(mapping: &[u32]) -> bool {
    mapping
        .iter()
        .enumerate()
        .all(|(position, target)| *target as usize == position)
}

/// Neighbours within the fragment and planarity from the whole component.
struct Guard {
    neighbours: Vec<Vec<usize>>,
    planar: Vec<bool>,
}

impl Guard {
    fn new(component: &Component, fragment: &Component) -> Self {
        let position = |name: &str| {
            fragment
                .atoms
                .iter()
                .position(|atom| atom.name.as_ref() == name)
        };
        let mut neighbours = vec![Vec::new(); fragment.atoms.len()];
        for bond in fragment.bonds.iter() {
            if let (Some(first), Some(second)) = (position(&bond.atom_a), position(&bond.atom_b)) {
                neighbours[first].push(second);
                neighbours[second].push(first);
            }
        }
        let planar = fragment
            .atoms
            .iter()
            .map(|atom| {
                component.bonds.iter().any(|bond| {
                    (bond.atom_a == atom.name || bond.atom_b == atom.name)
                        && (bond.aromatic
                            || !matches!(bond.order, BondOrder::Single | BondOrder::Unknown))
                })
            })
            .collect();
        Self { neighbours, planar }
    }

    fn admits(&self, mapping: &[u32]) -> bool {
        mapping.iter().enumerate().all(|(atom, target)| {
            *target as usize == atom
                || self.neighbours[atom]
                    .iter()
                    .any(|centre| mapping[*centre] as usize == *centre && self.planar[*centre])
        })
    }
}

#[cfg(test)]
#[path = "symmetry_tests.rs"]
mod tests;
