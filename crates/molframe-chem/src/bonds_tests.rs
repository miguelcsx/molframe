//! The grid search must find exactly the pairs a pairwise comparison finds.

use super::*;
use molframe_bench::{Sample, structure};
use std::collections::BTreeMap;

type Edge = (u32, u32);

/// Whether the pair may bond at this squared distance, asked of the two atoms alone.
fn candidate_allowed(
    left: AtomRef<'_>,
    right: AtomRef<'_>,
    thresholds: &[Option<f32>],
    distance_squared: f32,
) -> bool {
    let (Some(left_threshold), Some(right_threshold)) = (
        thresholds.get(left.index().as_usize()).copied().flatten(),
        thresholds.get(right.index().as_usize()).copied().flatten(),
    ) else {
        return false;
    };
    in_bonding_window(left_threshold, right_threshold, distance_squared)
        && pair_compatible(left, right)
}

/// Every bond the output holds, with its order and provenance as text.
fn edges(structure: &Structure) -> BTreeMap<Edge, String> {
    structure
        .data()
        .bonds
        .iter()
        .map(|bond| {
            (
                (bond.atom_a.get(), bond.atom_b.get()),
                format!("{:?}/{:?}", bond.order, bond.provenance),
            )
        })
        .collect()
}

/// The bonds a pairwise search over every atom pair would add.
fn pairwise(structure: &Structure) -> BTreeMap<Edge, String> {
    let mut output = BondTableBuilder::new();
    let mut existing = HashSet::new();
    for bond in structure.data().bonds.iter() {
        existing.insert(endpoints(bond.atom_a.get(), bond.atom_b.get()));
        output.push(bond);
    }
    add_polymer_links(structure, &mut output, &mut existing);
    let thresholds: Vec<_> = structure
        .data()
        .atoms()
        .map(|atom| atom.element().and_then(element_threshold))
        .collect();
    let positions = structure.positions();
    let count = u32::try_from(positions.len()).expect("a fixture fits in u32");
    for first in 0..count {
        for second in first + 1..count {
            let (Some(left), Some(right)) = (
                structure.atom(AtomIndex::new(first)),
                structure.atom(AtomIndex::new(second)),
            ) else {
                continue;
            };
            let distance = squared_distance(positions[first as usize], positions[second as usize]);
            if candidate_allowed(left, right, &thresholds, distance)
                && existing.insert(endpoints(first, second))
            {
                output.push(BondRecord {
                    atom_a: left.index(),
                    atom_b: right.index(),
                    order: intra_residue_order(left, right),
                    provenance: BondProvenance::InferredDistance,
                });
            }
        }
    }
    let mut data = structure.data().clone();
    data.bonds = output.finish();
    edges(&Structure::new(data))
}

#[test]
fn the_grid_finds_the_bonds_a_pairwise_search_finds() {
    for sample in [Sample::Tiny, Sample::Small, Sample::Medium] {
        let structure = structure(sample);
        let found = perceive_bonds(&structure).expect("coordinates are dense");
        let expected = pairwise(&structure);
        assert!(!expected.is_empty(), "{sample:?} has bonds");
        assert_eq!(edges(&found), expected, "{sample:?}");
    }
}

#[test]
fn perception_is_deterministic_and_keeps_file_bonds_first() {
    let structure = structure(Sample::Small);
    let first = perceive_bonds(&structure).expect("coordinates are dense");
    let second = perceive_bonds(&structure).expect("coordinates are dense");
    assert_eq!(edges(&first), edges(&second));
    let before = edges(&structure);
    let after = edges(&first);
    for (edge, description) in before {
        assert_eq!(after.get(&edge), Some(&description), "{edge:?}");
    }
}

/// Every pair of atoms in neighbouring cells, found through `grid`.
fn neighbouring_pairs(grid: &CellGrid<GridAtom>) -> Vec<Edge> {
    let mut pairs = Vec::new();
    grid.for_each_cell(|own, neighbourhood| {
        for first in &grid.items()[own] {
            for range in neighbourhood {
                for second in &grid.items()[range.clone()] {
                    if second.atom > first.atom {
                        pairs.push((first.atom, second.atom));
                    }
                }
            }
        }
    });
    pairs.sort_unstable();
    pairs
}

fn entries_of(
    positions: &[[f32; 3]],
    thresholds: &[Option<f32>],
    cutoff: f32,
) -> Vec<(crate::grid::CellKey, GridAtom)> {
    let mut entries = Vec::new();
    for (atom, (position, threshold)) in positions.iter().zip(thresholds).enumerate() {
        if let (Some(threshold), Some(cell)) = (*threshold, cell_for(*position, cutoff)) {
            entries.push((
                cell,
                GridAtom {
                    atom: u32::try_from(atom).expect("fits"),
                    position: *position,
                    threshold,
                },
            ));
        }
    }
    entries
}

#[test]
fn the_dense_and_sparse_grids_visit_the_same_neighbouring_pairs() {
    let mut state = 0x9e37_79b9_7f4a_7c15_u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        f32::from(u16::try_from(state >> 48).expect("16 bits fit")) / 65_536.0
    };
    let positions: Vec<[f32; 3]> = (0..1500)
        .map(|_| [next() * 60.0 - 30.0, next() * 45.0, next() * 60.0 - 10.0])
        .collect();
    let thresholds: Vec<Option<f32>> = (0..positions.len())
        .map(|index| (index % 11 != 0).then_some(1.75))
        .collect();
    let entries = entries_of(&positions, &thresholds, 2.77);
    let dense = CellGrid::dense(&entries).expect("a 60 A box is small");
    let sparse = CellGrid::sparse(entries.clone());
    let from_dense = neighbouring_pairs(&dense);
    assert_eq!(from_dense, neighbouring_pairs(&sparse));
    let mut brute = Vec::new();
    for (index, (left_cell, left)) in entries.iter().enumerate() {
        for (right_cell, right) in &entries[index + 1..] {
            let near = (0..3).all(|axis| (left_cell.0[axis] - right_cell.0[axis]).abs() <= 1);
            if near {
                brute.push((left.atom.min(right.atom), left.atom.max(right.atom)));
            }
        }
    }
    brute.sort_unstable();
    assert_eq!(from_dense, brute);
}

#[test]
fn a_span_too_large_to_tabulate_falls_back_to_the_sparse_grid() {
    let positions = [
        [0.0, 0.0, 0.0],
        [1.5, 0.0, 0.0],
        [90_000.0, 90_000.0, 90_000.0],
    ];
    let thresholds = [Some(1.75), Some(1.75), Some(1.75)];
    let grid = CellGrid::build(entries_of(&positions, &thresholds, 2.77));
    assert!(!grid.is_dense());
    assert_eq!(neighbouring_pairs(&grid), vec![(0, 1)]);
}
