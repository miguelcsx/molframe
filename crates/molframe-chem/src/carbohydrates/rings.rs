//! Bounded, conformer-consistent five/six-membered O/C cycle perception.

use super::geometry::ring_geometry;
use super::snfg::UNKNOWN;
use super::types::{Monosaccharide, SnfgSymbol};
use crate::Component;
use molframe_core::structure::ResidueRef;
use molframe_core::{AtomIndex, BondAdjacency, Element, Structure};

pub(super) fn residue_rings(
    structure: &Structure,
    graph: &BondAdjacency,
    residue: ResidueRef<'_>,
    component: Option<&Component>,
    symbol: Option<SnfgSymbol>,
    output: &mut Vec<Monosaccharide>,
) {
    let search = RingSearch {
        structure,
        graph,
        residue,
        component,
        symbol,
    };
    let mut path = [AtomIndex::new(0); 6];
    for atom in residue
        .atoms()
        .filter(|atom| atom.element() == Some(Element::OXYGEN))
    {
        path[0] = atom.index();
        search.extend(&mut path, 1, output);
    }
}

struct RingSearch<'a> {
    structure: &'a Structure,
    graph: &'a BondAdjacency,
    residue: ResidueRef<'a>,
    component: Option<&'a Component>,
    symbol: Option<SnfgSymbol>,
}

impl RingSearch<'_> {
    fn extend(&self, path: &mut [AtomIndex; 6], length: usize, output: &mut Vec<Monosaccharide>) {
        let current = path[length - 1];
        // The caller bounds malformed high-degree graphs before traversal.
        for next in self.graph.neighbours(current) {
            if *next == path[0] {
                if length >= 5 && path[1] < path[length - 1] && self.chordless(&path[..length]) {
                    self.emit(&path[..length], output);
                }
                continue;
            }
            if length == 6 || path[..length].contains(next) {
                continue;
            }
            let Some(atom) = self.structure.atom(*next) else {
                continue;
            };
            if atom.element() != Some(Element::CARBON)
                || atom
                    .residue()
                    .map(molframe_core::structure::ResidueRef::index)
                    != Some(self.residue.index())
                || !path[..length]
                    .iter()
                    .all(|previous| compatible(self.structure, *previous, *next))
            {
                continue;
            }
            path[length] = *next;
            self.extend(path, length + 1, output);
        }
    }

    fn chordless(&self, ring: &[AtomIndex]) -> bool {
        for (i, atom) in ring.iter().enumerate() {
            for (j, other) in ring.iter().enumerate().skip(i + 1) {
                if j == i + 1 || (i == 0 && j == ring.len() - 1) {
                    continue;
                }
                if self.graph.neighbours(*atom).contains(other) {
                    return false;
                }
            }
        }
        true
    }

    fn emit(&self, ring: &[AtomIndex], output: &mut Vec<Monosaccharide>) {
        let mut ring = ring.to_vec();
        let first = self.is_anomeric(ring[1], &ring);
        let last = self.is_anomeric(ring[ring.len() - 1], &ring);
        let anomeric_atom = match (first, last) {
            (true, false) => Some(ring[1]),
            (false, true) => {
                ring[1..].reverse();
                Some(ring[1])
            }
            _ => None,
        };
        output.push(Monosaccharide {
            residue: self.residue.index(),
            geometry: ring_geometry(self.structure, &ring, anomeric_atom),
            ring_atoms: ring,
            anomeric_atom,
            symbol: match self.symbol {
                Some(symbol) => symbol,
                None => UNKNOWN,
            },
        });
    }

    fn is_anomeric(&self, carbon: AtomIndex, ring: &[AtomIndex]) -> bool {
        let Some(atom) = self.structure.atom(carbon) else {
            return false;
        };
        if let (Some(component), Some(name)) = (self.component, atom.name()) {
            // The dictionary still identifies an anomeric site when its leaving
            // hydroxyl is absent from a glycosylated experimental structure.
            for bond in component.bonds.iter() {
                let other = if bond.atom_a.as_ref() == name {
                    &bond.atom_b
                } else if bond.atom_b.as_ref() == name {
                    &bond.atom_a
                } else {
                    continue;
                };
                let Some(expected) = component.atom(other) else {
                    continue;
                };
                if hetero(expected.element)
                    && !ring.iter().any(|index| {
                        self.structure
                            .atom(*index)
                            .and_then(molframe_core::structure::AtomRef::name)
                            == Some(other.as_ref())
                    })
                {
                    return true;
                }
            }
        }
        self.graph.neighbours(carbon).iter().any(|other| {
            !ring.contains(other)
                && self
                    .structure
                    .atom(*other)
                    .and_then(molframe_core::structure::AtomRef::element)
                    .is_some_and(hetero)
        })
    }
}

pub(super) fn compatible(structure: &Structure, a: AtomIndex, b: AtomIndex) -> bool {
    match (
        structure
            .atom(a)
            .and_then(molframe_core::structure::AtomRef::alt_id),
        structure
            .atom(b)
            .and_then(molframe_core::structure::AtomRef::alt_id),
    ) {
        (Some(a), Some(b)) => a.is_blank() || b.is_blank() || a == b,
        _ => true,
    }
}

pub(super) const fn hetero(element: Element) -> bool {
    matches!(
        element,
        Element::OXYGEN | Element::NITROGEN | Element::SULFUR
    )
}
