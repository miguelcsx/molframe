use crate::cif_fixture::read_cif;
use std::fmt::Write;

use super::{HeavyAtomHydrogenBondOptions, HeavyAtomHydrogenBondTable, heavy_atom_hydrogen_bonds};
use crate::hbond::{HeavyRoleTable, HydrogenBondError};
use molframe_bench::{Sample, structure};
use molframe_chem::{TableIdentity, TableKind};
use molframe_core::{
    AnnotationColumn, AtomAnnotation, AtomIndex, ExecutionContext, Presence, Structure,
};
use molframe_spatial::SpatialBackend;

const HEADER: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n\
_atom_site.label_asym_id\n_atom_site.label_seq_id\n\
_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";

/// One atom row: element, atom name, altloc (`.` for none), residue name,
/// residue number and coordinates.
type Row = (
    &'static str,
    &'static str,
    &'static str,
    &'static str,
    u32,
    [f32; 3],
);

fn build(rows: &[Row]) -> Structure {
    let mut source = String::from(HEADER);
    for (serial, (element, name, alt, residue, number, [x, y, z])) in rows.iter().enumerate() {
        let _ = writeln!(
            source,
            "ATOM {} {element} {name} {alt} {residue} A {number} {x} {y} {z}",
            serial + 1
        );
    }
    read_cif(&source)
}

/// The role table shipped under `data/`, decoded once; the repository carries it.
fn roles() -> &'static HeavyRoleTable {
    static ROLES: std::sync::OnceLock<HeavyRoleTable> = std::sync::OnceLock::new();
    ROLES.get_or_init(|| {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../data/heavy-atom-roles.roles"
        );
        let Ok(bytes) = std::fs::read(path) else {
            panic!("the repository's {path} is missing");
        };
        match HeavyRoleTable::decode(&bytes, TableIdentity::supplied("", "", "", "")) {
            Ok(provided) => provided.table,
            Err(error) => panic!("{path}: {error}"),
        }
    })
}

fn run(structure: &Structure, options: HeavyAtomHydrogenBondOptions) -> HeavyAtomHydrogenBondTable {
    match heavy_atom_hydrogen_bonds(structure, roles(), options, &ExecutionContext::default()) {
        Ok(table) => table,
        Err(error) => panic!("valid request failed: {error}"),
    }
}

fn pairs(table: &HeavyAtomHydrogenBondTable) -> Vec<(u32, u32)> {
    table
        .donor()
        .iter()
        .zip(table.acceptor())
        .map(|(donor, acceptor)| (donor.get(), acceptor.get()))
        .collect()
}

/// Ser OG donating to the carbonyl O of the next residue at an idealised 2.8 A.
fn ser_to_carbonyl(oxygen: [f32; 3], carbon: [f32; 3]) -> Structure {
    build(&[
        ("O", "OG", ".", "SER", 1, [0.0, 0.0, 0.0]),
        ("O", "O", ".", "GLY", 2, oxygen),
        ("C", "C", ".", "GLY", 2, carbon),
    ])
}

#[test]
fn a_serine_hydroxyl_to_a_carbonyl_oxygen_at_2_8_angstrom_is_one_directed_bond() {
    let structure = ser_to_carbonyl([2.8, 0.0, 0.0], [3.4, 1.07, 0.0]);
    let table = run(&structure, HeavyAtomHydrogenBondOptions::default());
    assert_eq!(pairs(&table), vec![(0, 1)]);
    assert!((table.distance()[0] - 2.8).abs() < 1.0e-5);
    // N...O=C: cos = -0.6 / |(0.6, 1.07)|.
    let expected = (-0.6_f64 / (0.6_f64.powi(2) + 1.07_f64.powi(2)).sqrt())
        .acos()
        .to_degrees();
    assert!((table.acceptor_angle_degrees()[0] - expected).abs() < 1.0e-3);
}

#[test]
fn distance_window_is_inclusive_of_the_cutoff_and_rejects_clashes() {
    let far = ser_to_carbonyl([3.6, 0.0, 0.0], [4.2, 1.07, 0.0]);
    assert!(run(&far, HeavyAtomHydrogenBondOptions::default()).is_empty());
    let wide = HeavyAtomHydrogenBondOptions {
        maximum_distance: 3.7,
        ..HeavyAtomHydrogenBondOptions::default()
    };
    assert_eq!(run(&far, wide).len(), 1);
    let close = ser_to_carbonyl([2.0, 0.0, 0.0], [2.6, 1.07, 0.0]);
    assert!(run(&close, HeavyAtomHydrogenBondOptions::default()).is_empty());
    let permissive = HeavyAtomHydrogenBondOptions {
        minimum_distance: 1.5,
        ..HeavyAtomHydrogenBondOptions::default()
    };
    assert_eq!(run(&close, permissive).len(), 1);
}

#[test]
fn an_acute_donor_acceptor_antecedent_angle_is_rejected_unless_unfiltered() {
    // The carbon sits on the donor side of the oxygen: the donor approaches
    // from behind the carbonyl bond, about 49 degrees.
    let structure = ser_to_carbonyl([2.8, 0.0, 0.0], [2.2, 0.7, 0.0]);
    assert!(run(&structure, HeavyAtomHydrogenBondOptions::default()).is_empty());
    let unfiltered = HeavyAtomHydrogenBondOptions {
        minimum_acceptor_angle_degrees: None,
        ..HeavyAtomHydrogenBondOptions::default()
    };
    let table = run(&structure, unfiltered);
    assert_eq!(table.len(), 1);
    assert!(table.acceptor_angle_degrees()[0] < 90.0);
}

#[test]
fn an_acceptor_without_antecedents_reports_nan_and_is_kept() {
    let structure = build(&[
        ("O", "O", ".", "HOH", 1, [0.0, 0.0, 0.0]),
        ("O", "O", ".", "HOH", 2, [2.8, 0.0, 0.0]),
    ]);
    let table = run(&structure, HeavyAtomHydrogenBondOptions::default());
    // Waters donate and accept, so both orientations are reported.
    assert_eq!(pairs(&table), vec![(0, 1), (1, 0)]);
    assert!(table.acceptor_angle_degrees().iter().all(|a| a.is_nan()));
}

#[test]
fn a_peptide_neighbour_pair_is_removed_by_covalent_adjacency_not_by_name() {
    // O(i)=C(i)-N(i+1): a 1-3 pair about 2.24 A apart.
    let structure = build(&[
        ("O", "O", ".", "GLY", 1, [0.0, 0.0, 0.0]),
        ("C", "C", ".", "GLY", 1, [1.23, 0.0, 0.0]),
        ("N", "N", ".", "GLY", 2, [1.93, 1.13, 0.0]),
    ]);
    let relaxed = HeavyAtomHydrogenBondOptions {
        minimum_distance: 0.0,
        minimum_acceptor_angle_degrees: None,
        ..HeavyAtomHydrogenBondOptions::default()
    };
    assert!(run(&structure, relaxed).is_empty());
    // With no atoms regarded as bonded the very same pair is accepted.
    let no_bonds = HeavyAtomHydrogenBondOptions {
        covalent_distance: 0.5,
        ..relaxed
    };
    assert_eq!(pairs(&run(&structure, no_bonds)), vec![(2, 0)]);
}

#[test]
fn pairs_inside_one_residue_are_dropped_unless_requested() {
    let structure = build(&[
        ("O", "OG", ".", "SER", 1, [0.0, 0.0, 0.0]),
        ("O", "O", ".", "SER", 1, [2.8, 0.0, 0.0]),
    ]);
    let options = HeavyAtomHydrogenBondOptions::default();
    assert!(run(&structure, options).is_empty());
    let same = HeavyAtomHydrogenBondOptions {
        include_same_residue_pairs: true,
        ..options
    };
    assert_eq!(pairs(&run(&structure, same)), vec![(0, 1)]);
}

#[test]
fn alternate_locations_pair_only_when_compatible() {
    let count = |donor, acceptor| {
        let structure = build(&[
            ("O", "OG", donor, "SER", 1, [0.0, 0.0, 0.0]),
            ("O", "O", acceptor, "GLY", 2, [2.8, 0.0, 0.0]),
        ]);
        run(&structure, HeavyAtomHydrogenBondOptions::default()).len()
    };
    assert_eq!(count("A", "B"), 0);
    assert_eq!(count("A", "A"), 1);
    assert_eq!(count(".", "B"), 1);
    assert_eq!(count("A", "."), 1);
}

#[test]
fn proline_amide_nitrogen_does_not_donate_except_at_the_chain_start() {
    let rows = |first: &'static str| {
        build(&[
            ("N", "N", ".", first, 1, [0.0, 0.0, 0.0]),
            ("O", "O", ".", "HOH", 2, [2.9, 0.0, 0.0]),
        ])
    };
    let options = HeavyAtomHydrogenBondOptions::default();
    // The water accepts, so a donating N yields one bond; otherwise only
    // the water as donor to nothing.
    assert_eq!(pairs(&run(&rows("PRO"), options)), vec![(0, 1)]);
    let internal = build(&[
        ("C", "C", ".", "GLY", 1, [-3.0, 3.0, 0.0]),
        ("N", "N", ".", "PRO", 2, [0.0, 0.0, 0.0]),
        ("O", "O", ".", "HOH", 3, [2.9, 0.0, 0.0]),
    ]);
    assert!(run(&internal, options).is_empty());
}

#[test]
fn sulfur_is_typed_only_on_request() {
    let structure = build(&[
        ("S", "SG", ".", "CYS", 1, [0.0, 0.0, 0.0]),
        ("O", "O", ".", "HOH", 2, [3.3, 0.0, 0.0]),
    ]);
    assert!(run(&structure, HeavyAtomHydrogenBondOptions::default()).is_empty());
    let with_sulfur = HeavyAtomHydrogenBondOptions {
        include_sulfur: true,
        ..HeavyAtomHydrogenBondOptions::default()
    };
    assert_eq!(run(&structure, with_sulfur).len(), 2);
}

#[test]
fn side_chain_roles_follow_the_standard_tables() {
    // Lys NZ (donor only) against Asp OD1 (acceptor only): one direction.
    let structure = build(&[
        ("N", "NZ", ".", "LYS", 1, [0.0, 0.0, 0.0]),
        ("O", "OD1", ".", "ASP", 2, [2.8, 0.0, 0.0]),
    ]);
    assert_eq!(
        pairs(&run(&structure, HeavyAtomHydrogenBondOptions::default())),
        vec![(0, 1)]
    );
    // Two acceptors never bond.
    let acceptors = build(&[
        ("O", "OD1", ".", "ASP", 1, [0.0, 0.0, 0.0]),
        ("O", "OE1", ".", "GLU", 2, [2.8, 0.0, 0.0]),
    ]);
    assert!(run(&acceptors, HeavyAtomHydrogenBondOptions::default()).is_empty());
    // Nucleic acids: guanine N1 donates to cytosine N3.
    let pair = build(&[
        ("N", "N1", ".", "DG", 1, [0.0, 0.0, 0.0]),
        ("N", "N3", ".", "DC", 2, [2.9, 0.0, 0.0]),
    ]);
    assert_eq!(
        pairs(&run(&pair, HeavyAtomHydrogenBondOptions::default())),
        vec![(0, 1)]
    );
}

fn with_roles(structure: &Structure, donors: &[u32], acceptors: &[u32]) -> Structure {
    let column = |set: &[u32]| {
        AnnotationColumn::from_entries((0..structure.atom_count()).map(|atom| {
            if set.contains(&atom) {
                (true, Presence::Present)
            } else {
                (false, Presence::Inapplicable)
            }
        }))
        .expect("small annotation column")
    };
    let mut data = structure.data().clone();
    data.annotations.insert(
        molframe_core::HBOND_DONOR_ANNOTATION,
        AtomAnnotation::Boolean(column(donors)),
    );
    data.annotations.insert(
        molframe_core::HBOND_ACCEPTOR_ANNOTATION,
        AtomAnnotation::Boolean(column(acceptors)),
    );
    Structure::new(data)
}

#[test]
fn unknown_groups_use_n_and_o_both_ways_until_annotations_say_otherwise() {
    let structure = build(&[
        ("N", "N1", ".", "LIG", 1, [0.0, 0.0, 0.0]),
        ("O", "O1", ".", "LIG", 2, [2.9, 0.0, 0.0]),
        ("C", "C1", ".", "LIG", 3, [0.0, 2.9, 0.0]),
    ]);
    let options = HeavyAtomHydrogenBondOptions::default();
    assert_eq!(pairs(&run(&structure, options)), vec![(0, 1), (1, 0)]);
    let annotated = with_roles(&structure, &[0], &[1]);
    assert_eq!(pairs(&run(&annotated, options)), vec![(0, 1)]);
}

#[test]
fn invalid_options_are_refused() {
    let structure = build(&[("O", "O", ".", "HOH", 1, [0.0, 0.0, 0.0])]);
    let bad = [
        HeavyAtomHydrogenBondOptions {
            maximum_distance: f32::NAN,
            ..HeavyAtomHydrogenBondOptions::default()
        },
        HeavyAtomHydrogenBondOptions {
            minimum_distance: 4.0,
            ..HeavyAtomHydrogenBondOptions::default()
        },
        HeavyAtomHydrogenBondOptions {
            minimum_acceptor_angle_degrees: Some(200.0),
            ..HeavyAtomHydrogenBondOptions::default()
        },
        HeavyAtomHydrogenBondOptions {
            covalent_distance: 0.0,
            ..HeavyAtomHydrogenBondOptions::default()
        },
    ];
    for options in bad {
        assert_eq!(
            heavy_atom_hydrogen_bonds(&structure, roles(), options, &ExecutionContext::default())
                .err(),
            Some(HydrogenBondError::InvalidOptions)
        );
    }
}

fn atom_name(structure: &Structure, atom: AtomIndex) -> String {
    match structure
        .data()
        .atom(atom)
        .and_then(molframe_core::structure::AtomRef::name)
    {
        Some(name) => name.to_owned(),
        None => panic!("atom without a name"),
    }
}

fn residue_number(structure: &Structure, atom: AtomIndex) -> i64 {
    let Some(residue) = structure
        .data()
        .atom(atom)
        .and_then(molframe_core::structure::AtomRef::residue)
    else {
        panic!("atom without a residue");
    };
    i64::from(residue.index().get())
}

#[test]
fn crambin_main_chain_hydrogen_bonds_match_the_literature_range() {
    let crambin = structure(Sample::Tiny);
    let table = run(&crambin, HeavyAtomHydrogenBondOptions::default());
    let backbone: Vec<_> = table
        .donor()
        .iter()
        .zip(table.acceptor())
        .filter(|&(&donor, &acceptor)| {
            atom_name(&crambin, donor) == "N" && atom_name(&crambin, acceptor) == "O"
        })
        .map(|(&donor, &acceptor)| {
            (
                residue_number(&crambin, donor),
                residue_number(&crambin, acceptor),
            )
        })
        .collect();
    // Crambin (46 residues, two helices and a short sheet) has roughly 30 to 45
    // main-chain N-H...O bonds under HBPLUS-style criteria.
    assert!(
        (30..=45).contains(&backbone.len()),
        "main-chain bonds: {}",
        backbone.len()
    );
    // Helix i -> i+4: the first helix (residues 7-19) contributes NH(i)...O(i-4).
    let helical = backbone.iter().filter(|&&(n, o)| n - o == 4).count();
    assert!(helical >= 8, "i,i+4 bonds: {helical}");
    println!(
        "crambin backbone {} helical {helical} total {}",
        backbone.len(),
        table.len()
    );
}

#[test]
fn results_are_identical_across_backends_and_worker_counts() {
    let haemoglobin = structure(Sample::Medium);
    let reference = run(&haemoglobin, HeavyAtomHydrogenBondOptions::default());
    assert!(!reference.is_empty());
    for backend in [
        SpatialBackend::BruteForce,
        SpatialBackend::CellList,
        SpatialBackend::KdTree,
    ] {
        for workers in [1, 4] {
            let Ok(context) = ExecutionContext::builder().worker_budget(workers).build() else {
                panic!("context");
            };
            let options = HeavyAtomHydrogenBondOptions {
                backend,
                ..HeavyAtomHydrogenBondOptions::default()
            };
            let Ok(table) = heavy_atom_hydrogen_bonds(&haemoglobin, roles(), options, &context)
            else {
                panic!("request failed");
            };
            assert_eq!(pairs(&table), pairs(&reference), "{backend:?} x {workers}");
        }
    }
    let sorted = pairs(&reference);
    assert!(sorted.windows(2).all(|w| w[0] < w[1]));
    // Regression golden over haemoglobin (four chains, hem and waters) with
    // the default policy; it guards typing and exclusions, not a reference tool.
    assert_eq!(reference.len(), 1501);
}

/// Roles of the two atoms of a pair: the number of directed bonds each
/// orientation yields, `(first donates, second donates)`.
fn directions(first: Row, second: Row, options: HeavyAtomHydrogenBondOptions) -> (bool, bool) {
    let found = pairs(&run(&build(&[first, second]), options));
    (found.contains(&(0, 1)), found.contains(&(1, 0)))
}

#[test]
fn histidine_tautomers_donate_from_the_protonated_nitrogen_only() {
    let water = ("O", "O", ".", "HOH", 2, [2.9, 0.0, 0.0]);
    let options = HeavyAtomHydrogenBondOptions::default();
    let nitrogen = |residue: &'static str, name: &'static str| {
        directions(
            ("N", name, ".", residue, 1, [0.0, 0.0, 0.0]),
            water,
            options,
        )
    };
    // Water does both: a donor-only nitrogen takes no bond from it and an
    // acceptor-only nitrogen gives none to it.
    assert_eq!(nitrogen("HID", "ND1"), (true, false));
    assert_eq!(nitrogen("HID", "NE2"), (false, true));
    assert_eq!(nitrogen("HIE", "NE2"), (true, false));
    assert_eq!(nitrogen("HIE", "ND1"), (false, true));
    assert_eq!(nitrogen("HIP", "ND1"), (true, false));
    assert_eq!(nitrogen("HIP", "NE2"), (true, false));
}

#[test]
fn neutral_acids_and_amines_gain_the_roles_their_state_implies() {
    let options = HeavyAtomHydrogenBondOptions::default();
    let donor_probe = |residue: &'static str, atom: &'static str, element: &'static str| {
        // Against Asp OD1 (acceptor only): can the variant atom donate?
        directions(
            (element, atom, ".", residue, 1, [0.0, 0.0, 0.0]),
            ("O", "OD1", ".", "ASP", 2, [2.8, 0.0, 0.0]),
            options,
        )
        .0
    };
    assert!(donor_probe("LYS", "NZ", "N"));
    assert!(donor_probe("LYN", "NZ", "N"));
    assert!(donor_probe("ASH", "OD2", "O"));
    assert!(!donor_probe("ASH", "OD1", "O"));
    assert!(donor_probe("GLH", "OE2", "O"));
    assert!(!donor_probe("GLH", "OE1", "O"));
    assert!(!donor_probe("ASP", "OD2", "O"));
}

#[test]
fn disulphide_and_thiolate_sulfur_accept_but_never_donate() {
    let with_sulfur = HeavyAtomHydrogenBondOptions {
        include_sulfur: true,
        ..HeavyAtomHydrogenBondOptions::default()
    };
    for residue in ["CYX", "CYM"] {
        let sulfur = ("S", "SG", ".", residue, 1, [0.0, 0.0, 0.0]);
        let donor = ("N", "NZ", ".", "LYS", 2, [3.3, 0.0, 0.0]);
        assert_eq!(directions(sulfur, donor, with_sulfur), (false, true));
        assert!(
            run(
                &build(&[sulfur, donor]),
                HeavyAtomHydrogenBondOptions::default()
            )
            .is_empty()
        );
    }
}
