use crate::smarts::{AtomTest, BondExpression, SmartsDataError, SmartsMatch, SmartsPattern};
use crate::{Component, StereoConfiguration};
use pdbiox_core::annotation::{
    AROMATIC_ATOM_ANNOTATION, AtomAnnotation, FORMAL_CHARGE_ANNOTATION,
    STEREO_CONFIGURATION_ANNOTATION,
};
use pdbiox_core::topology::{Csr, CsrBuilder};
use pdbiox_core::{BondOrder, Element, Structure};
use std::collections::{BTreeMap, VecDeque};

mod matcher;

use matcher::Matcher;

pub(crate) fn find_matches(pattern: &SmartsPattern, component: &Component) -> Vec<SmartsMatch> {
    let graph = Graph::from_component(component);
    Matcher::new(pattern, &graph).find(None)
}

pub(crate) fn has_match(pattern: &SmartsPattern, component: &Component) -> bool {
    let graph = Graph::from_component(component);
    Matcher::new(pattern, &graph).has(None)
}

pub(crate) fn find_structure_matches(
    pattern: &SmartsPattern,
    structure: &Structure,
) -> Result<Vec<SmartsMatch>, SmartsDataError> {
    validate_structure_data(pattern, structure)?;
    let graph = Graph::from_structure(structure)?;
    Ok(Matcher::new(pattern, &graph).find(None))
}

#[derive(Clone, Copy)]
struct GraphAtom {
    element: Element,
    charge: Option<i8>,
    aromatic: Option<bool>,
    stereo: Option<StereoConfiguration>,
}

#[derive(Clone, Copy)]
struct GraphBond {
    order: BondOrder,
}

struct Graph {
    atoms: Vec<GraphAtom>,
    bonds: Vec<GraphBond>,
    adjacency: Csr<(usize, usize)>,
}

impl Graph {
    fn from_component(component: &Component) -> Self {
        let atoms = component
            .atoms
            .iter()
            .map(|atom| GraphAtom {
                element: atom.element,
                charge: Some(atom.charge),
                aromatic: Some(atom.aromatic),
                stereo: atom.stereo,
            })
            .collect();
        let names: BTreeMap<&str, usize> = component
            .atoms
            .iter()
            .enumerate()
            .map(|(index, atom)| (atom.name.as_ref(), index))
            .collect();
        let edges = component.bonds.iter().filter_map(|bond| {
            let first = *names.get(bond.atom_a.as_ref())?;
            let second = *names.get(bond.atom_b.as_ref())?;
            Some((first, second, GraphBond { order: bond.order }))
        });
        Self::build(atoms, edges)
    }

    fn from_structure(structure: &Structure) -> Result<Self, SmartsDataError> {
        if !structure.data().bonds.is_available() {
            return Err(SmartsDataError::Connectivity);
        }
        let data = structure.data();
        let aromatic = data.annotations.get(AROMATIC_ATOM_ANNOTATION);
        let charges = data.annotations.get(FORMAL_CHARGE_ANNOTATION);
        let stereo = data.annotations.get(STEREO_CONFIGURATION_ANNOTATION);
        let atoms = data
            .atoms()
            .map(|atom| {
                let index = atom.index().get();
                GraphAtom {
                    element: match atom.element() {
                        Some(element) => element,
                        None => Element::UNKNOWN,
                    },
                    charge: integer_annotation(charges, index)
                        .and_then(|value| i8::try_from(value).ok()),
                    aromatic: boolean_annotation(aromatic, index),
                    stereo: symbol_annotation(stereo, index, structure).and_then(stereo_value),
                }
            })
            .collect();
        let edges = data.bonds.iter().map(|bond| {
            (
                bond.atom_a.as_usize(),
                bond.atom_b.as_usize(),
                GraphBond { order: bond.order },
            )
        });
        Ok(Self::build(atoms, edges))
    }

    fn build(
        atoms: Vec<GraphAtom>,
        edges: impl IntoIterator<Item = (usize, usize, GraphBond)>,
    ) -> Self {
        // Edges are collected once so their degrees can size a flat adjacency,
        // which replaces one heap allocation per atom with two shared buffers.
        let mut bonds = Vec::new();
        let mut accepted = Vec::new();
        let mut degrees = vec![0usize; atoms.len()];
        for (first, second, bond) in edges {
            if first >= atoms.len() || second >= atoms.len() {
                continue;
            }
            let edge = bonds.len();
            bonds.push(bond);
            accepted.push((first, second, edge));
            degrees[first] += 1;
            degrees[second] += 1;
        }

        let mut builder = CsrBuilder::with_degrees(&degrees);
        for (first, second, edge) in accepted {
            builder.push(first, (second, edge));
            builder.push(second, (first, edge));
        }
        let mut adjacency = builder.finish();
        for atom in 0..adjacency.rows() {
            adjacency.row_mut(atom).sort_unstable();
        }

        Self {
            atoms,
            bonds,
            adjacency,
        }
    }

    fn bond(&self, first: usize, second: usize) -> Option<GraphBond> {
        let neighbours = self.adjacency.row(first);
        let upper = neighbours.partition_point(|(other, _)| *other <= second);
        let (other, edge) = neighbours.get(upper.checked_sub(1)?)?;
        if *other != second {
            return None;
        }
        self.bonds.get(*edge).copied()
    }

    fn edge_in_ring(&self, first: usize, second: usize) -> bool {
        self.path_avoiding(first, second, usize::MAX, ordered(first, second))
            .is_some()
    }

    fn ring_bond_count(&self, atom: usize) -> usize {
        self.adjacency
            .row(atom)
            .iter()
            .filter(|(other, _)| self.edge_in_ring(atom, *other))
            .count()
    }

    fn smallest_ring(&self, atom: usize) -> Option<usize> {
        let neighbours = self.adjacency.row(atom);
        let mut smallest = None;
        for left in 0..neighbours.len() {
            for right in left + 1..neighbours.len() {
                let excluded = ordered(atom, neighbours[left].0);
                let Some(path) =
                    self.path_avoiding(neighbours[left].0, neighbours[right].0, atom, excluded)
                else {
                    continue;
                };
                smallest = Some(smallest.map_or(path + 2, |value: usize| value.min(path + 2)));
            }
        }
        smallest
    }

    fn ring_count(&self, atom: usize) -> usize {
        let bonds = self.ring_bond_count(atom);
        if bonds == 0 { 0 } else { bonds - 1 }
    }

    fn path_avoiding(
        &self,
        start: usize,
        goal: usize,
        avoided: usize,
        excluded: (usize, usize),
    ) -> Option<usize> {
        let mut distances = vec![usize::MAX; self.adjacency.rows()];
        let mut queue = VecDeque::from([start]);
        distances[start] = 0;
        while let Some(atom) = queue.pop_front() {
            if atom == goal {
                return Some(distances[atom]);
            }
            for (next, _) in self.adjacency.row(atom) {
                if *next == avoided
                    || ordered(atom, *next) == excluded
                    || distances[*next] != usize::MAX
                {
                    continue;
                }
                distances[*next] = distances[atom] + 1;
                queue.push_back(*next);
            }
        }
        None
    }
}

fn validate_structure_data(
    pattern: &SmartsPattern,
    structure: &Structure,
) -> Result<(), SmartsDataError> {
    if !structure.data().bonds.is_available() {
        return Err(SmartsDataError::Connectivity);
    }
    let annotations = &structure.data().annotations;
    if pattern_uses(pattern, |test| {
        matches!(test, AtomTest::Aromatic | AtomTest::Aliphatic)
    }) && annotations.get(AROMATIC_ATOM_ANNOTATION).is_none()
    {
        return Err(SmartsDataError::Aromaticity);
    }
    if pattern_uses(pattern, |test| matches!(test, AtomTest::Charge(_)))
        && annotations.get(FORMAL_CHARGE_ANNOTATION).is_none()
    {
        return Err(SmartsDataError::FormalCharge);
    }
    if pattern_uses(pattern, |test| matches!(test, AtomTest::Stereo(_)))
        && annotations.get(STEREO_CONFIGURATION_ANNOTATION).is_none()
    {
        return Err(SmartsDataError::Stereochemistry);
    }
    Ok(())
}

fn pattern_uses(pattern: &SmartsPattern, predicate: impl Fn(&AtomTest) -> bool + Copy) -> bool {
    pattern.atoms.iter().any(|expression| {
        expression.alternatives.iter().flatten().any(|signed| {
            predicate(&signed.test)
                || matches!(&signed.test, AtomTest::Recursive(nested) if pattern_uses(nested, predicate))
        })
    })
}

fn boolean_annotation(annotation: Option<&AtomAnnotation>, atom: u32) -> Option<bool> {
    let AtomAnnotation::Boolean(column) = annotation? else {
        return None;
    };
    let (value, presence) = column.get(atom)?;
    presence.is_present().then_some(value)
}

fn integer_annotation(annotation: Option<&AtomAnnotation>, atom: u32) -> Option<i64> {
    let AtomAnnotation::Integer(column) = annotation? else {
        return None;
    };
    let (value, presence) = column.get(atom)?;
    presence.is_present().then_some(value)
}

fn symbol_annotation<'a>(
    annotation: Option<&AtomAnnotation>,
    atom: u32,
    structure: &'a Structure,
) -> Option<&'a str> {
    let AtomAnnotation::Symbol(column) = annotation? else {
        return None;
    };
    let (symbol, presence) = column.get(atom)?;
    presence
        .is_present()
        .then(|| structure.resolve(symbol))
        .flatten()
}

fn stereo_value(value: &str) -> Option<StereoConfiguration> {
    match value {
        "R" => Some(StereoConfiguration::R),
        "S" => Some(StereoConfiguration::S),
        "mixed" => Some(StereoConfiguration::Mixed),
        _ => None,
    }
}

fn bond_matches(expression: BondExpression, bond: GraphBond, in_ring: bool) -> bool {
    match expression {
        BondExpression::Default => matches!(
            bond.order,
            BondOrder::Single | BondOrder::Unknown | BondOrder::Aromatic
        ),
        BondExpression::Any => true,
        BondExpression::Order(order) => bond.order == order,
        BondExpression::Ring(expected) => in_ring == expected,
        BondExpression::NotOrder(order) => bond.order != order,
    }
}

fn valence(graph: &Graph, atom: usize) -> u16 {
    graph
        .adjacency
        .row(atom)
        .iter()
        .filter_map(|(other, _)| graph.bond(atom, *other))
        .map(|bond| match bond.order {
            BondOrder::Single | BondOrder::Polymeric | BondOrder::Unknown | BondOrder::Aromatic => {
                1
            }
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Quadruple => 4,
        })
        .sum()
}

fn ordered(first: usize, second: usize) -> (usize, usize) {
    if first < second {
        (first, second)
    } else {
        (second, first)
    }
}

#[cfg(test)]
#[path = "smarts_tests.rs"]
mod tests;
