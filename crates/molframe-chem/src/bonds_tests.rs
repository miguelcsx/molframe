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
    for bond in structure.data().bonds.iter() {
        output.push(bond);
    }
    add_polymer_links(structure, &mut output);
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
            if candidate_allowed(left, right, &thresholds, distance) {
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

#[test]
fn chunked_cell_visits_reproduce_the_whole_sequence() {
    let mut positions: Vec<[f32; 3]> = (0..400)
        .map(|n| {
            let n = f32::from(u16::try_from(n).expect("small"));
            [n.rem_euclid(17.0) * 1.3, n.rem_euclid(11.0) * 1.7, n * 0.05]
        })
        .collect();
    let thresholds = vec![Some(1.75); positions.len()];
    let dense_entries = entries_of(&positions, &thresholds, 2.77);
    positions.push([90_000.0, 90_000.0, 90_000.0]);
    let sparse_entries = entries_of(&positions, &vec![Some(1.75); positions.len()], 2.77);
    for entries in [dense_entries, sparse_entries] {
        let grid = CellGrid::build(entries);
        let mut whole = Vec::new();
        grid.for_each_cell(|own, neighbourhood| whole.push((own, neighbourhood.to_vec())));
        for block in [1, 7, 64, grid.cell_count().max(1)] {
            let mut pieces = Vec::new();
            let mut start = 0;
            while start < grid.cell_count() {
                let end = (start + block).min(grid.cell_count());
                grid.visit_cells(start..end, |own, neighbourhood| {
                    pieces.push((own, neighbourhood.to_vec()));
                });
                start = end;
            }
            assert_eq!(pieces, whole, "block size {block}");
        }
    }
}

#[test]
fn perception_is_identical_at_every_worker_budget() {
    let structure = structure(Sample::Large);
    let mut reference = None;
    for workers in [1, 2, 4, 16] {
        let context = ExecutionContext::builder()
            .worker_budget(workers)
            .build()
            .expect("a positive worker budget is valid");
        let found = perceive_bonds_in(&structure, &context).expect("coordinates are dense");
        let observed = edges(&found);
        match &reference {
            None => reference = Some(observed),
            Some(first) => assert_eq!(&observed, first, "{workers} workers"),
        }
    }
}

fn atom_pair(elements: [&str; 2], alternates: [&str; 2], distance: f32) -> Structure {
    let source = format!(
        "data_pair\nloop_\n_atom_site.id\n_atom_site.type_symbol\n\
         _atom_site.label_atom_id\n_atom_site.label_comp_id\n\
         _atom_site.label_asym_id\n_atom_site.label_seq_id\n\
         _atom_site.label_alt_id\n_atom_site.Cartn_x\n\
         _atom_site.Cartn_y\n_atom_site.Cartn_z\n\
         1 {} QL LIG A 1 {} 0 0 0\n2 {} QR LIG A 1 {} {distance} 0 0\n",
        elements[0], alternates[0], elements[1], alternates[1],
    );
    let input = molframe_core::InputBuffer::from_bytes(source.into_bytes());
    molframe_cif::read(&input, &molframe_core::ReadOptions::new())
        .expect("a two-atom fixture parses")
        .0
}

#[test]
fn the_public_covalent_window_is_inclusive_and_rejects_nonfinite_distances() {
    let source = atom_pair(["C", "C"], [".", "."], 1.0);
    let left = source.atom(AtomIndex::new(0)).expect("left atom");
    let right = source.atom(AtomIndex::new(1)).expect("right atom");
    let minimum = 0.16_f32;
    let maximum_distance = 3.50_f32 / 1.95;
    let maximum = maximum_distance * maximum_distance;
    for distance in [minimum, minimum.next_up(), maximum.next_down(), maximum] {
        assert!(covalent_pair(left, right, distance), "{distance}");
        assert!(covalent_pair(right, left, distance), "symmetric {distance}");
    }
    for distance in [
        minimum.next_down(),
        maximum.next_up(),
        -1.0,
        f32::NAN,
        f32::INFINITY,
        f32::NEG_INFINITY,
    ] {
        assert!(!covalent_pair(left, right, distance), "{distance}");
    }
}

#[test]
fn perception_and_public_covalent_eligibility_share_chemical_boundaries() {
    let cases = [
        (["C", "O"], [".", "."], 1.2, true, true),
        (["C", "C"], [".", "."], 0.399, false, false),
        (["C", "C"], [".", "."], 0.4, true, true),
        (["C", "C"], [".", "."], 1.79, true, true),
        (["C", "C"], [".", "."], 1.80, false, false),
        (["H", "H"], [".", "."], 0.7, false, false),
        (["H", "C"], [".", "."], 1.0, true, true),
        (["C", "O"], ["A", "B"], 1.2, false, false),
        (["C", "O"], ["A", "A"], 1.2, true, true),
        (["C", "O"], [".", "B"], 1.2, true, true),
        (["B", "C"], [".", "."], 1.2, false, false),
        (["?", "C"], [".", "."], 1.2, false, false),
        (["Fe", "N"], [".", "."], 2.0, false, true),
        (["Fe", "N"], ["A", "B"], 2.0, false, false),
        (["Fe", "N"], [".", "."], 2.3, false, false),
    ];
    for (elements, alternates, distance, covalent, perceived) in cases {
        let source = atom_pair(elements, alternates, distance);
        let left = source.atom(AtomIndex::new(0)).expect("left atom");
        let right = source.atom(AtomIndex::new(1)).expect("right atom");
        let label = format!("{elements:?} {alternates:?} {distance}");
        assert_eq!(
            covalent_pair(left, right, distance * distance),
            covalent,
            "{label}"
        );
        assert_eq!(
            covalent_pair(right, left, distance * distance),
            covalent,
            "{label}"
        );
        let found = perceive_bonds(&source).expect("dense coordinates");
        let bonds: Vec<_> = found.data().bonds.iter().collect();
        assert_eq!(bonds.len(), usize::from(perceived), "{label}");
        if perceived {
            assert_eq!(bonds[0].atom_a, left.index(), "{label}");
            assert_eq!(bonds[0].atom_b, right.index(), "{label}");
            assert_eq!(bonds[0].order, BondOrder::Single, "{label}");
            assert_eq!(
                bonds[0].provenance,
                BondProvenance::InferredDistance,
                "{label}"
            );
        }
    }
}
