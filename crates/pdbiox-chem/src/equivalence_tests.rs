use super::*;
use crate::{ComponentAtom, ComponentKind};
use pdbiox_core::{BondOrder, Element};
use std::collections::BTreeMap;

fn atom(name: &str, element: Element) -> ComponentAtom {
    ComponentAtom {
        name: name.into(),
        alternate_name: None,
        element,
        charge: 0,
        aromatic: false,
        leaving: false,
        stereo: None,
    }
}

fn bond(atom_a: &str, atom_b: &str, order: BondOrder) -> ComponentBond {
    ComponentBond {
        atom_a: atom_a.into(),
        atom_b: atom_b.into(),
        order,
        aromatic: order == BondOrder::Aromatic,
        stereo: None,
    }
}

fn component() -> Component {
    Component {
        id: "ASP".into(),
        name: "ASPARTATE".into(),
        kind: ComponentKind::AminoAcid,
        parent: None,
        one_letter_code: Some(b'D'),
        formula: None,
        atoms: vec![
            atom("CG", Element::CARBON),
            atom("OD1", Element::OXYGEN),
            atom("OD2", Element::OXYGEN),
            atom("CB", Element::CARBON),
        ]
        .into(),
        bonds: vec![
            bond("CG", "OD1", BondOrder::Aromatic),
            bond("CG", "OD2", BondOrder::Aromatic),
            bond("CG", "CB", BondOrder::Single),
        ]
        .into(),
        ideal_coordinates: None,
        model_coordinates: None,
    }
}

#[test]
fn exact_automorphisms_are_complete_ordered_and_bounded() {
    let component = component();
    let mappings = automorphisms(&component, 4)
        .unwrap_or_else(|error| panic!("automorphism enumeration failed: {error}"));
    assert_eq!(mappings.len(), 2);
    assert_eq!(mappings[0].as_ref(), &[0, 1, 2, 3]);
    assert_eq!(mappings[1].as_ref(), &[0, 2, 1, 3]);
    assert_eq!(
        automorphisms(&component, 1),
        Err(AutomorphismLimit { limit: 1 })
    );
}

#[test]
fn equivalent_atoms_are_derived_from_graph_symmetry_not_names() {
    let classes = equivalence_classes(&component());
    assert!(classes.equivalent(1, 2));
    assert!(!classes.equivalent(0, 3));
    assert_eq!(
        classes.classes(),
        &[Arc::from([0]), Arc::from([1, 2]), Arc::from([3])]
    );
}

#[test]
fn changing_one_bond_order_breaks_the_oxygen_equivalence() {
    let mut component = component();
    Arc::make_mut(&mut component.bonds)[0].order = BondOrder::Double;
    Arc::make_mut(&mut component.bonds)[0].aromatic = false;
    let classes = equivalence_classes(&component);
    assert!(!classes.equivalent(1, 2));
}

#[test]
fn cache_is_explicit_versioned_and_reuses_the_computed_result() {
    let mut cache = EquivalenceCache::new(DictionaryVersion::new("2026-08-01"));
    let component = component();
    let first = cache.get(&component);
    let second = cache.get(&component);
    assert!(Arc::ptr_eq(&first, &second));
    assert_eq!(cache.len(), 1);
    assert_eq!(cache.version().as_str(), "2026-08-01");
}

#[test]
fn refined_backtracking_matches_exhaustive_permutations_on_all_four_atom_graphs() {
    const COUNT: usize = 4;
    for element_mask in 0u16..1 << COUNT {
        let elements: Vec<u8> = (0..COUNT)
            .map(|atom| {
                if element_mask & (1 << atom) == 0 {
                    6
                } else {
                    8
                }
            })
            .collect();
        for edge_mask in 0u16..1 << 6 {
            let mut edges = BTreeMap::new();
            let mut bit = 0;
            for atom_a in 0..COUNT {
                for atom_b in atom_a + 1..COUNT {
                    if edge_mask & (1 << bit) != 0 {
                        edges.insert((atom_a, atom_b), 1);
                        edges.insert((atom_b, atom_a), 1);
                    }
                    bit += 1;
                }
            }
            let graph = Graph {
                colours: refine_colours(&elements, &edges),
                elements: elements.clone(),
                edges,
            };
            for source in 0..COUNT {
                for target in 0..COUNT {
                    assert_eq!(
                        graph.has_automorphism(source, target),
                        brute_automorphism(&graph, source, target),
                        "element_mask={element_mask}, edge_mask={edge_mask}, {source}->{target}"
                    );
                }
            }
        }
    }
}

fn brute_automorphism(graph: &Graph, source: usize, target: usize) -> bool {
    let mut permutation: Vec<usize> = (0..graph.elements.len()).collect();
    permutations(&mut permutation, 0, &mut |candidate| {
        candidate[source] == target
            && candidate.iter().enumerate().all(|(atom, mapped)| {
                graph.elements[atom] == graph.elements[*mapped]
                    && (0..candidate.len()).all(|other| {
                        graph.edge(atom, other) == graph.edge(*mapped, candidate[other])
                    })
            })
    })
}

fn permutations(
    values: &mut [usize],
    position: usize,
    predicate: &mut impl FnMut(&[usize]) -> bool,
) -> bool {
    if position == values.len() {
        return predicate(values);
    }
    for swap in position..values.len() {
        values.swap(position, swap);
        if permutations(values, position + 1, predicate) {
            return true;
        }
        values.swap(position, swap);
    }
    false
}
