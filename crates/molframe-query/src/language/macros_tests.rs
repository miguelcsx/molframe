use super::*;
use crate::{Groups, Query, col};
use molframe_core::contract::AnalysisPolicy;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::{AtomIndex, BondOrder, BondProvenance, BondRecord, BondTableBuilder};

const SOURCE: &str = "data_hydrogens
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 N N GLY A 1 0 0 0
ATOM 2 O O GLY A 1 10 0 0
ATOM 3 S S GLY A 1 20 0 0
ATOM 4 C C GLY A 1 30 0 0
ATOM 5 P P GLY A 1 40 0 0
ATOM 6 F F GLY A 1 50 0 0
ATOM 7 H HN GLY A 1 60 0 0
ATOM 8 H HO GLY A 1 70 0 0
ATOM 9 H HS GLY A 1 80 0 0
ATOM 10 H HC GLY A 1 90 0 0
ATOM 11 H FREE GLY A 1 100 0 0
ATOM 12 H MIXED GLY A 1 110 0 0
ATOM 13 H HP GLY A 1 120 0 0
ATOM 14 H HF GLY A 1 130 0 0
ATOM 15 H HH1 GLY A 1 140 0 0
ATOM 16 H HH2 GLY A 1 150 0 0
";

fn raw_structure() -> Structure {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    molframe_cif::read(&input, &ReadOptions::new())
        .expect("fixture parses without inferred bonds")
        .0
}

fn with_edges(structure: &Structure, edges: &[(u32, u32)]) -> Structure {
    let mut bonds = BondTableBuilder::new();
    for &(a, b) in edges {
        bonds.push(BondRecord {
            atom_a: AtomIndex::new(a),
            atom_b: AtomIndex::new(b),
            order: BondOrder::Single,
            provenance: BondProvenance::File,
        });
    }
    let mut data = structure.data().clone();
    data.bonds = bonds.finish();
    Structure::new(data)
}

fn structure() -> Structure {
    with_edges(
        &raw_structure(),
        &[
            (0, 6),
            (1, 7),
            (2, 8),
            (3, 9),
            (3, 11),
            (0, 11),
            (4, 12),
            (5, 13),
            (14, 15),
        ],
    )
}

fn selected(structure: &Structure, source: &str, universe: &AtomSelection) -> AtomSelection {
    Query::compile(source)
        .expect("query compiles")
        .evaluate_in(
            structure,
            universe,
            &AnalysisPolicy::default(),
            &Groups::new(),
            None,
        )
        .expect("bond topology is available")
        .selection
}

#[test]
fn hydrogen_polarity_is_an_exhaustive_partition_by_actual_bonded_elements() {
    let structure = structure();
    let all = AtomSelection::All(structure.atom_count());
    let polar = selected(&structure, "polar_hydrogen", &all);
    let nonpolar = selected(&structure, "nonpolar_hydrogen", &all);
    assert_eq!(polar.iter().collect::<Vec<_>>(), [6, 7, 8, 11]);
    assert_eq!(nonpolar.iter().collect::<Vec<_>>(), [9, 10, 12, 13, 14, 15]);
    assert!(polar.intersect(&nonpolar).is_empty());
    assert_eq!(
        polar.union(&nonpolar),
        selected(&structure, "hydrogen", &all)
    );
    assert_eq!(
        selected(&structure, "hydrogen and not polar_hydrogen", &all),
        nonpolar
    );
}

#[test]
fn hydrogen_polarity_reads_parents_outside_restricted_universes_and_conjunctions() {
    let structure = structure();
    let universe = AtomSelection::from_sorted(vec![6, 7, 9, 10]);
    assert_eq!(
        selected(&structure, "polar_hydrogen", &universe)
            .iter()
            .collect::<Vec<_>>(),
        [6, 7]
    );
    assert_eq!(
        selected(&structure, "nonpolar_hydrogen", &universe)
            .iter()
            .collect::<Vec<_>>(),
        [9, 10]
    );
    for source in ["hydrogen and polar_hydrogen", "polar_hydrogen and hydrogen"] {
        assert_eq!(
            selected(&structure, source, &universe)
                .iter()
                .collect::<Vec<_>>(),
            [6, 7]
        );
    }
    let all = AtomSelection::All(structure.atom_count());
    assert_eq!(
        selected(&structure, "index 6 and polar_hydrogen", &all)
            .iter()
            .collect::<Vec<_>>(),
        [6]
    );
}

#[test]
fn hydrogen_polarity_follows_topology_not_atom_names_or_distances() {
    let raw = raw_structure();
    let nitrogen_parent = with_edges(&raw, &[(0, 9)]);
    let carbon_parent = with_edges(&raw, &[(3, 9)]);
    let universe = AtomSelection::from_sorted(vec![9]);
    assert_eq!(
        selected(&nitrogen_parent, "polar_hydrogen", &universe),
        universe
    );
    assert!(selected(&carbon_parent, "polar_hydrogen", &universe).is_empty());
    assert_eq!(
        selected(&carbon_parent, "nonpolar_hydrogen", &universe),
        universe
    );
}

#[test]
fn known_empty_bonds_make_all_explicit_hydrogen_nonpolar() {
    let structure = with_edges(&raw_structure(), &[]);
    let all = AtomSelection::All(structure.atom_count());
    assert!(selected(&structure, "polar_hydrogen", &all).is_empty());
    assert_eq!(
        selected(&structure, "nonpolar_hydrogen", &all),
        selected(&structure, "hydrogen", &all)
    );
}

#[test]
fn unavailable_bonds_error_instead_of_claiming_hydrogens_are_unbonded() {
    let structure = raw_structure();
    assert!(!structure.data().bonds.is_available());
    for name in ["polar_hydrogen", "nonpolar_hydrogen"] {
        let findings = Query::compile(name)
            .expect("query compiles")
            .evaluate(&structure, &AnalysisPolicy::default(), &Groups::new(), None)
            .expect_err("unavailable connectivity must fail");
        assert_eq!(findings[0].code(), Code::E4003);
        let error = macro_selection(
            &structure,
            &AtomSelection::Empty,
            Macro::from_name(name).expect("macro exists"),
            &mut Vec::new(),
        )
        .expect_err("empty output does not establish topology availability");
        assert_eq!(error.code(), Code::E4003);
    }
}

#[test]
fn hydrogen_builders_text_and_roundtrips_share_the_same_query_and_result() {
    let structure = structure();
    for builder in [col::polar_hydrogen(), col::nonpolar_hydrogen()] {
        let query = Query::from_builder(builder);
        let printed = crate::print::print(&query.logical_plan().0);
        let parsed = Query::compile(&printed).expect("printed query parses");
        assert_eq!(query.fingerprint(), parsed.fingerprint());
        let evaluate = |query: &Query| {
            query
                .evaluate(&structure, &AnalysisPolicy::default(), &Groups::new(), None)
                .expect("topology is available")
                .selection
        };
        assert_eq!(evaluate(&query), evaluate(&parsed));
    }
}
