use crate::smarts::{AtomTest, BondExpression, SmartsDataError, SmartsMatch, SmartsPattern};
use crate::{Component, StereoConfiguration};
use pdbiox_core::annotation::{
    AROMATIC_ATOM_ANNOTATION, AtomAnnotation, FORMAL_CHARGE_ANNOTATION,
    STEREO_CONFIGURATION_ANNOTATION,
};
use pdbiox_core::{BondOrder, Element, Structure};
use std::collections::{BTreeMap, BTreeSet, VecDeque};

pub(crate) fn find_matches(pattern: &SmartsPattern, component: &Component) -> Vec<SmartsMatch> {
    let graph = Graph::from_component(component);
    Matcher::new(pattern, &graph).find(None)
}

pub(crate) fn has_match(pattern: &SmartsPattern, component: &Component) -> bool {
    let graph = Graph::from_component(component);
    !Matcher::new(pattern, &graph).find(None).is_empty()
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
    adjacency: Vec<Vec<(usize, usize)>>,
    edge_by_pair: BTreeMap<(usize, usize), usize>,
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
        let mut bonds = Vec::new();
        let mut adjacency = vec![Vec::new(); atoms.len()];
        let mut edge_by_pair = BTreeMap::new();
        for (first, second, bond) in edges {
            if first >= atoms.len() || second >= atoms.len() {
                continue;
            }
            let edge = bonds.len();
            bonds.push(bond);
            adjacency[first].push((second, edge));
            adjacency[second].push((first, edge));
            edge_by_pair.insert(ordered(first, second), edge);
        }
        for neighbours in &mut adjacency {
            neighbours.sort_unstable();
        }
        Self {
            atoms,
            bonds,
            adjacency,
            edge_by_pair,
        }
    }

    fn bond(&self, first: usize, second: usize) -> Option<GraphBond> {
        self.edge_by_pair
            .get(&ordered(first, second))
            .and_then(|edge| self.bonds.get(*edge))
            .copied()
    }

    fn edge_in_ring(&self, first: usize, second: usize) -> bool {
        self.path_avoiding(first, second, usize::MAX, ordered(first, second))
            .is_some()
    }

    fn ring_bond_count(&self, atom: usize) -> usize {
        self.adjacency[atom]
            .iter()
            .filter(|(other, _)| self.edge_in_ring(atom, *other))
            .count()
    }

    fn smallest_ring(&self, atom: usize) -> Option<usize> {
        let neighbours = &self.adjacency[atom];
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
        let mut distances = vec![usize::MAX; self.adjacency.len()];
        let mut queue = VecDeque::from([start]);
        distances[start] = 0;
        while let Some(atom) = queue.pop_front() {
            if atom == goal {
                return Some(distances[atom]);
            }
            for (next, _) in &self.adjacency[atom] {
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

struct Matcher<'a> {
    pattern: &'a SmartsPattern,
    graph: &'a Graph,
    pattern_adjacency: Vec<Vec<(usize, BondExpression)>>,
}

impl<'a> Matcher<'a> {
    fn new(pattern: &'a SmartsPattern, graph: &'a Graph) -> Self {
        let mut pattern_adjacency = vec![Vec::new(); pattern.atoms.len()];
        for bond in &pattern.bonds {
            pattern_adjacency[bond.first].push((bond.second, bond.expression));
            pattern_adjacency[bond.second].push((bond.first, bond.expression));
        }
        Self {
            pattern,
            graph,
            pattern_adjacency,
        }
    }

    fn find(&self, anchored_first: Option<usize>) -> Vec<SmartsMatch> {
        if self.pattern.atoms.is_empty() {
            return Vec::new();
        }
        let mut mappings = BTreeSet::new();
        let candidates: Box<dyn Iterator<Item = usize>> = match anchored_first {
            Some(atom) => Box::new(std::iter::once(atom)),
            None => Box::new(0..self.graph.atoms.len()),
        };
        for atom in candidates {
            if !self.atom_matches(0, atom) {
                continue;
            }
            let mut mapping = vec![None; self.pattern.atoms.len()];
            let mut used = BTreeSet::from([atom]);
            mapping[0] = Some(atom);
            self.search(&mut mapping, &mut used, &mut mappings);
        }
        mappings
            .into_iter()
            .map(|indices| SmartsMatch {
                atom_indices: indices.into_boxed_slice(),
            })
            .collect()
    }

    fn search(
        &self,
        mapping: &mut [Option<usize>],
        used: &mut BTreeSet<usize>,
        results: &mut BTreeSet<Vec<usize>>,
    ) {
        let Some(query) = self.next_query_atom(mapping) else {
            results.insert(mapping.iter().flatten().copied().collect());
            return;
        };
        for candidate in 0..self.graph.atoms.len() {
            if used.contains(&candidate)
                || !self.atom_matches(query, candidate)
                || !self.connected_constraints_match(query, candidate, mapping)
            {
                continue;
            }
            mapping[query] = Some(candidate);
            used.insert(candidate);
            self.search(mapping, used, results);
            used.remove(&candidate);
            mapping[query] = None;
        }
    }

    fn next_query_atom(&self, mapping: &[Option<usize>]) -> Option<usize> {
        (0..mapping.len())
            .filter(|query| mapping[*query].is_none())
            .max_by_key(|query| {
                self.pattern_adjacency[*query]
                    .iter()
                    .filter(|(other, _)| mapping[*other].is_some())
                    .count()
            })
    }

    fn connected_constraints_match(
        &self,
        query: usize,
        candidate: usize,
        mapping: &[Option<usize>],
    ) -> bool {
        self.pattern_adjacency[query]
            .iter()
            .all(|(other, expression)| {
                let Some(mapped_other) = mapping[*other] else {
                    return true;
                };
                self.graph
                    .bond(candidate, mapped_other)
                    .is_some_and(|bond| {
                        bond_matches(
                            *expression,
                            bond,
                            self.graph.edge_in_ring(candidate, mapped_other),
                        )
                    })
            })
    }

    fn atom_matches(&self, query: usize, atom: usize) -> bool {
        self.pattern.atoms[query]
            .alternatives
            .iter()
            .any(|alternative| {
                alternative
                    .iter()
                    .all(|signed| self.test_matches(&signed.test, atom) != signed.negated)
            })
    }

    fn test_matches(&self, test: &AtomTest, atom: usize) -> bool {
        let graph_atom = self.graph.atoms[atom];
        match test {
            AtomTest::Any => true,
            AtomTest::Element(element) => graph_atom.element == *element,
            AtomTest::Aromatic => graph_atom.aromatic == Some(true),
            AtomTest::Aliphatic => graph_atom.aromatic == Some(false),
            AtomTest::Degree(degree) | AtomTest::Connectivity(degree) => {
                self.graph.adjacency[atom].len() == usize::from(*degree)
            }
            AtomTest::Hydrogens(count) => self.hydrogen_count(atom) == usize::from(*count),
            AtomTest::Charge(charge) => graph_atom.charge == Some(*charge),
            AtomTest::RingCount(None) => self.graph.ring_count(atom) > 0,
            AtomTest::RingCount(Some(count)) => self.graph.ring_count(atom) == usize::from(*count),
            AtomTest::RingSize(None) => self.graph.smallest_ring(atom).is_some(),
            AtomTest::RingSize(Some(size)) => {
                self.graph.smallest_ring(atom) == Some(usize::from(*size))
            }
            AtomTest::Valence(value) => valence(self.graph, atom) == u16::from(*value),
            AtomTest::RingBonds(count) => self.graph.ring_bond_count(atom) == usize::from(*count),
            AtomTest::Stereo(stereo) => graph_atom.stereo == Some(*stereo),
            AtomTest::Recursive(pattern) => !Matcher::new(pattern, self.graph)
                .find(Some(atom))
                .is_empty(),
        }
    }

    fn hydrogen_count(&self, atom: usize) -> usize {
        self.graph.adjacency[atom]
            .iter()
            .filter(|(other, _)| self.graph.atoms[*other].element == Element::HYDROGEN)
            .count()
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
    graph.adjacency[atom]
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
