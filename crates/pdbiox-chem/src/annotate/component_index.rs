//! A component's name index and its precomputed hydrogen-bonding roles.
//!
//! `Component::atom` is a linear scan comparing atom names, and the
//! hydrogen-bonding predicates call it from inside nested bond loops: deciding
//! whether one nitrogen is an amide walks every bond, and for each of those
//! walks every bond again, resolving a name linearly each time. That is
//! `O(atoms² × bonds²)` string comparisons for a single observed atom.
//!
//! Both results depend only on the component and the atom name, never on the
//! residue, so they are identical for every glycine in a structure. Computing
//! them once per component and indexing by name turns the per-atom cost into a
//! binary search.

use crate::{Component, ComponentAtom};
use pdbiox_core::{BondOrder, Element};
use std::sync::Arc;

/// A component with its name index and per-atom hydrogen-bonding roles.
#[derive(Debug)]
pub(crate) struct ComponentIndex {
    /// The component every position below refers to.
    component: Arc<Component>,
    /// Atom positions ordered by name, so a name is found by binary search.
    order: Box<[u32]>,
    /// Whether the atom at each position can donate a hydrogen bond.
    donor: Box<[bool]>,
    /// Whether the atom at each position can accept one.
    acceptor: Box<[bool]>,
}

impl ComponentIndex {
    /// Indexes `component` and evaluates both roles for every one of its atoms.
    pub(crate) fn build(component: Arc<Component>) -> Self {
        let count = component.atoms.len();
        let mut order: Vec<u32> = (0..count)
            .filter_map(|position| u32::try_from(position).ok())
            .collect();
        order.sort_unstable_by(|left, right| {
            name_at(&component, *left).cmp(name_at(&component, *right))
        });
        let order = order.into_boxed_slice();

        let donor: Vec<bool> = (0..count)
            .map(|position| evaluate_donor(&component, &order, position))
            .collect();
        let acceptor: Vec<bool> = (0..count)
            .map(|position| evaluate_acceptor(&component, &order, position))
            .collect();

        Self {
            component,
            order,
            donor: donor.into_boxed_slice(),
            acceptor: acceptor.into_boxed_slice(),
        }
    }

    /// The indexed component.
    pub(crate) fn component(&self) -> &Component {
        &self.component
    }

    /// The position of the atom called `name`, in `O(log atoms)`.
    pub(crate) fn position(&self, name: &str) -> Option<usize> {
        position_in(&self.component, &self.order, name)
    }

    /// The atom called `name`, or `None` when the component has no such atom.
    pub(crate) fn atom(&self, name: &str) -> Option<&ComponentAtom> {
        self.component.atoms.get(self.position(name)?)
    }

    /// Whether the named atom can donate a hydrogen bond.
    pub(crate) fn is_donor(&self, name: &str) -> bool {
        match self.position(name).and_then(|slot| self.donor.get(slot)) {
            Some(donor) => *donor,
            None => false,
        }
    }

    /// Whether the named atom can accept a hydrogen bond.
    pub(crate) fn is_acceptor(&self, name: &str) -> bool {
        match self.position(name).and_then(|slot| self.acceptor.get(slot)) {
            Some(acceptor) => *acceptor,
            None => false,
        }
    }
}

/// The name of the atom at `position`, or the empty string when out of range.
fn name_at(component: &Component, position: u32) -> &str {
    match component.atoms.get(position as usize) {
        Some(atom) => atom.name.as_ref(),
        None => "",
    }
}

/// Binary searches the name-ordered positions.
fn position_in(component: &Component, order: &[u32], name: &str) -> Option<usize> {
    let slot = order
        .binary_search_by(|position| name_at(component, *position).cmp(name))
        .ok()?;
    Some(*order.get(slot)? as usize)
}

/// The positions bonded to `position`, resolved through the name index.
fn bonded_positions(component: &Component, order: &[u32], position: usize) -> Vec<usize> {
    let Some(atom) = component.atoms.get(position) else {
        return Vec::new();
    };
    let name = atom.name.as_ref();

    component
        .bonds
        .iter()
        .filter_map(|bond| {
            if bond.atom_a.as_ref() == name {
                position_in(component, order, &bond.atom_b)
            } else if bond.atom_b.as_ref() == name {
                position_in(component, order, &bond.atom_a)
            } else {
                None
            }
        })
        .collect()
}

/// Whether the atom at `position` can donate a hydrogen bond.
fn evaluate_donor(component: &Component, order: &[u32], position: usize) -> bool {
    let Some(atom) = component.atoms.get(position) else {
        return false;
    };
    if !matches!(
        atom.element,
        Element::NITROGEN | Element::OXYGEN | Element::SULFUR
    ) {
        return false;
    }

    bonded_positions(component, order, position)
        .into_iter()
        .filter_map(|neighbour| component.atoms.get(neighbour))
        .any(|neighbour| neighbour.element.is_hydrogen())
}

/// Whether the atom at `position` can accept a hydrogen bond.
fn evaluate_acceptor(component: &Component, order: &[u32], position: usize) -> bool {
    let Some(atom) = component.atoms.get(position) else {
        return false;
    };

    match atom.element {
        Element::OXYGEN | Element::SULFUR | Element::FLUORINE | Element::CHLORINE => {
            atom.charge <= 0
        }
        Element::NITROGEN => {
            atom.charge <= 0
                && !is_amide(component, order, position)
                && !(atom.aromatic
                    && bonded_positions(component, order, position)
                        .into_iter()
                        .filter_map(|neighbour| component.atoms.get(neighbour))
                        .any(|neighbour| neighbour.element.is_hydrogen()))
        }
        _ => false,
    }
}

/// Whether the nitrogen at `position` is bonded to a carbonyl or thiocarbonyl.
fn is_amide(component: &Component, order: &[u32], position: usize) -> bool {
    bonded_positions(component, order, position)
        .into_iter()
        .any(|neighbour| {
            let Some(carbon) = component.atoms.get(neighbour) else {
                return false;
            };
            if carbon.element != Element::CARBON {
                return false;
            }
            let carbon_name = carbon.name.as_ref();

            component.bonds.iter().any(|bond| {
                if bond.order != BondOrder::Double {
                    return false;
                }
                let other = if bond.atom_a.as_ref() == carbon_name {
                    position_in(component, order, &bond.atom_b)
                } else if bond.atom_b.as_ref() == carbon_name {
                    position_in(component, order, &bond.atom_a)
                } else {
                    None
                };
                other
                    .and_then(|slot| component.atoms.get(slot))
                    .is_some_and(|atom| matches!(atom.element, Element::OXYGEN | Element::SULFUR))
            })
        })
}

#[cfg(test)]
#[path = "component_index_tests.rs"]
mod tests;
