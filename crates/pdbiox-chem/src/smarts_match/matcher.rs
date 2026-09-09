use super::{
    AtomTest, BondExpression, Element, Graph, SmartsMatch, SmartsPattern, bond_matches, valence,
};
use std::ops::ControlFlow;

pub(super) struct Matcher<'a> {
    pattern: &'a SmartsPattern,
    graph: &'a Graph,
    pattern_adjacency: Vec<Vec<(usize, BondExpression)>>,
}

impl<'a> Matcher<'a> {
    pub(super) fn new(pattern: &'a SmartsPattern, graph: &'a Graph) -> Self {
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

    pub(super) fn find(&self, anchored_first: Option<usize>) -> Vec<SmartsMatch> {
        if self.pattern.atoms.is_empty() {
            return Vec::new();
        }
        let mut mappings = Vec::new();
        self.for_each_mapping(anchored_first, |mapping| {
            mappings.push(mapping.iter().flatten().copied().collect::<Vec<_>>());
            ControlFlow::Continue(())
        });
        mappings.sort_unstable();
        mappings.dedup();
        mappings
            .into_iter()
            .map(|indices| SmartsMatch {
                atom_indices: indices.into_boxed_slice(),
            })
            .collect()
    }

    pub(super) fn has(&self, anchored_first: Option<usize>) -> bool {
        if self.pattern.atoms.is_empty() {
            return false;
        }
        let mut found = false;
        self.for_each_mapping(anchored_first, |_| {
            found = true;
            ControlFlow::Break(())
        });
        found
    }

    pub(super) fn for_each_mapping(
        &self,
        anchored_first: Option<usize>,
        mut visitor: impl FnMut(&[Option<usize>]) -> ControlFlow<()>,
    ) {
        let mut mapping = vec![None; self.pattern.atoms.len()];
        let mut used = vec![false; self.graph.atoms.len()];
        if let Some(atom) = anchored_first {
            if atom < self.graph.atoms.len() {
                let _ = self.start_at(atom, &mut mapping, &mut used, &mut visitor);
            }
        } else {
            for atom in 0..self.graph.atoms.len() {
                if self
                    .start_at(atom, &mut mapping, &mut used, &mut visitor)
                    .is_break()
                {
                    return;
                }
            }
        }
    }

    fn start_at<V>(
        &self,
        atom: usize,
        mapping: &mut [Option<usize>],
        used: &mut [bool],
        visitor: &mut V,
    ) -> ControlFlow<()>
    where
        V: FnMut(&[Option<usize>]) -> ControlFlow<()>,
    {
        if !self.atom_matches(0, atom) {
            return ControlFlow::Continue(());
        }
        used[atom] = true;
        mapping[0] = Some(atom);
        let flow = self.search(mapping, used, visitor);
        used[atom] = false;
        mapping[0] = None;
        flow
    }

    fn search<V>(
        &self,
        mapping: &mut [Option<usize>],
        used: &mut [bool],
        visitor: &mut V,
    ) -> ControlFlow<()>
    where
        V: FnMut(&[Option<usize>]) -> ControlFlow<()>,
    {
        let Some(query) = self.next_query_atom(mapping) else {
            return visitor(mapping);
        };
        let mapped_neighbour = self.pattern_adjacency[query]
            .iter()
            .find_map(|(other, _)| mapping[*other]);
        if let Some(neighbour) = mapped_neighbour {
            for &(candidate, _) in self.graph.adjacency.row(neighbour) {
                if self
                    .try_candidate(query, candidate, mapping, used, visitor)
                    .is_break()
                {
                    return ControlFlow::Break(());
                }
            }
        } else {
            for candidate in 0..self.graph.atoms.len() {
                if self
                    .try_candidate(query, candidate, mapping, used, visitor)
                    .is_break()
                {
                    return ControlFlow::Break(());
                }
            }
        }
        ControlFlow::Continue(())
    }

    fn try_candidate<V>(
        &self,
        query: usize,
        candidate: usize,
        mapping: &mut [Option<usize>],
        used: &mut [bool],
        visitor: &mut V,
    ) -> ControlFlow<()>
    where
        V: FnMut(&[Option<usize>]) -> ControlFlow<()>,
    {
        if used[candidate]
            || !self.atom_matches(query, candidate)
            || !self.connected_constraints_match(query, candidate, mapping)
        {
            return ControlFlow::Continue(());
        }
        mapping[query] = Some(candidate);
        used[candidate] = true;
        let flow = self.search(mapping, used, visitor);
        used[candidate] = false;
        mapping[query] = None;
        flow
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
                        let in_ring = matches!(expression, BondExpression::Ring(_))
                            && self.graph.edge_in_ring(candidate, mapped_other);
                        bond_matches(*expression, bond, in_ring)
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
                self.graph.adjacency.row(atom).len() == usize::from(*degree)
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
            AtomTest::Recursive(pattern) => Matcher::new(pattern, self.graph).has(Some(atom)),
        }
    }

    fn hydrogen_count(&self, atom: usize) -> usize {
        self.graph
            .adjacency
            .row(atom)
            .iter()
            .filter(|(other, _)| self.graph.atoms[*other].element == Element::HYDROGEN)
            .count()
    }
}
