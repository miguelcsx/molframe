//! Conservative graph aromaticity: continuous p orbitals and 4n+2 electrons.
//!
//! Enumerates simple circuits up to 24 atoms, including fused-ring perimeters.
//! Neutral C/N with one double bond donate one electron; singly bonded neutral
//! N/O/S and anionic C donate two; cationic C donates zero. Exocyclic multiple
//! bonds, unsupported valences/elements and unknown charges are not inferred.
//! Explicit aromatic bonds remain authoritative. Geometry is not consulted.
//! Cost is O(V+E) plus circuit enumeration, bounded to 100,000 edge visits per
//! molecule. Unproved atoms remain unknown, including after either bound.

use molframe_core::topology::{Csr, CsrBuilder};
use molframe_core::{BondOrder, Element};

pub(crate) struct Aromaticity {
    pub(crate) atoms: Vec<bool>,
    pub(crate) bonds: Vec<bool>,
}

pub(crate) fn perceive(
    elements: &[Element],
    charges: &[Option<i8>],
    bonds: &[(usize, usize, BondOrder)],
) -> Aromaticity {
    let mut result = Aromaticity {
        atoms: vec![false; elements.len()],
        bonds: vec![false; bonds.len()],
    };
    let mut degrees = vec![0; elements.len()];
    for &(a, b, _) in bonds {
        degrees[a] += 1;
        degrees[b] += 1;
    }
    let mut builder = CsrBuilder::with_degrees(&degrees);
    for (edge, &(a, b, order)) in bonds.iter().enumerate() {
        builder.push(a, (b, edge));
        builder.push(b, (a, edge));
        if order == BondOrder::Aromatic {
            result.atoms[a] = true;
            result.atoms[b] = true;
            result.bonds[edge] = true;
        }
    }
    let adjacency = builder.finish();
    let mut search = Circuits {
        elements,
        charges,
        bonds,
        adjacency,
        result,
        path: Vec::with_capacity(24),
        edges: Vec::with_capacity(24),
        visited: vec![false; elements.len()],
        remaining: 100_000,
    };
    for start in 0..elements.len() {
        if search.remaining == 0 {
            break;
        }
        search.path.push(start);
        search.visited[start] = true;
        search.walk(start, start);
        search.visited[start] = false;
        search.path.clear();
    }
    search.result
}

struct Circuits<'a> {
    elements: &'a [Element],
    charges: &'a [Option<i8>],
    bonds: &'a [(usize, usize, BondOrder)],
    adjacency: Csr<(usize, usize)>,
    result: Aromaticity,
    path: Vec<usize>,
    edges: Vec<usize>,
    visited: Vec<bool>,
    remaining: usize,
}

impl Circuits<'_> {
    fn walk(&mut self, start: usize, atom: usize) {
        for offset in 0..self.adjacency.row(atom).len() {
            if self.remaining == 0 {
                return;
            }
            self.remaining -= 1;
            let (next, edge) = self.adjacency.row(atom)[offset];
            if !matches!(self.bonds[edge].2, BondOrder::Single | BondOrder::Double) {
                continue;
            }
            if next == start && self.path.len() >= 3 {
                self.edges.push(edge);
                if self.aromatic() {
                    for &site in &self.path {
                        self.result.atoms[site] = true;
                    }
                    for &bond in &self.edges {
                        self.result.bonds[bond] = true;
                    }
                }
                self.edges.pop();
            } else if next > start && !self.visited[next] && self.path.len() < 24 {
                self.visited[next] = true;
                self.path.push(next);
                self.edges.push(edge);
                self.walk(start, next);
                self.edges.pop();
                self.path.pop();
                self.visited[next] = false;
            }
        }
    }

    fn aromatic(&self) -> bool {
        let mut electrons = 0;
        for &atom in &self.path {
            let Some(charge) = self.charges[atom] else {
                return false;
            };
            let mut doubles = 0;
            let mut valence = 0;
            for &(_, edge) in self.adjacency.row(atom) {
                match self.bonds[edge].2 {
                    BondOrder::Single => valence += 1,
                    BondOrder::Double if self.edges.contains(&edge) => {
                        doubles += 1;
                        valence += 2;
                    }
                    _ => return false,
                }
            }
            let contribution = match (self.elements[atom], charge, doubles) {
                (Element::CARBON, 0, 1) | (Element::NITROGEN, 1, 1) if valence <= 4 => 1,
                (Element::NITROGEN, 0, 1) if valence <= 3 => 1,
                (Element::NITROGEN, 0, 0) | (Element::CARBON, -1, 0) if valence <= 3 => 2,
                (Element::OXYGEN | Element::SULFUR, 0, 0) if valence == 2 => 2,
                (Element::CARBON, 1, 0) if valence <= 3 => 0,
                _ => return false,
            };
            electrons += contribution;
        }
        electrons >= 2 && electrons % 4 == 2
    }
}

#[cfg(test)]
#[path = "aromaticity_tests.rs"]
mod tests;
