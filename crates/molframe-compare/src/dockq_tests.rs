use super::{DockQOptions, dockq};
use crate::CompareError;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;
use std::fmt::Write;

const HEADER: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";

fn structure(body: &str) -> Structure {
    let input = InputBuffer::from_bytes(format!("{HEADER}{body}").into_bytes());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

/// One residue: backbone N, CA, C, O at fixed offsets from `origin`.
struct Residue {
    chain: &'static str,
    number: u32,
    origin: [f64; 3],
    extra: Vec<(&'static str, &'static str, [f64; 3])>,
}

const BACKBONE: [(&str, &str, [f64; 3]); 4] = [
    ("N", "N", [0.0, 0.0, 0.0]),
    ("C", "CA", [1.5, 0.0, 0.0]),
    ("C", "C", [2.0, 1.4, 0.0]),
    ("O", "O", [3.0, 1.5, 1.0]),
];

fn residue(chain: &'static str, number: u32, origin: [f64; 3]) -> Residue {
    Residue {
        chain,
        number,
        origin,
        extra: Vec::new(),
    }
}

fn body(residues: &[Residue]) -> String {
    let mut text = String::new();
    let mut id = 0;
    for residue in residues {
        let atoms = BACKBONE.iter().chain(residue.extra.iter());
        for (element, name, offset) in atoms {
            id += 1;
            let _ = writeln!(
                text,
                "ATOM {id} {element} {name} ALA {} {} {} {} {}",
                residue.chain,
                residue.number,
                residue.origin[0] + offset[0],
                residue.origin[1] + offset[1],
                residue.origin[2] + offset[2],
            );
        }
    }
    text
}

fn options() -> DockQOptions {
    DockQOptions::published()
}

fn score(model: &[Residue], native: &[Residue], options: DockQOptions) -> super::DockQ {
    match dockq(
        &structure(&body(model)),
        &structure(&body(native)),
        "A",
        "B",
        options,
    ) {
        Ok(result) => result,
        Err(error) => panic!("valid score: {error}"),
    }
}

/// Two separate contact sites 20 Å apart; each ligand residue touches only the
/// receptor residue beneath it.
fn two_sites(second_ligand_height: f64) -> Vec<Residue> {
    vec![
        residue("A", 1, [0.0, 0.0, 0.0]),
        residue("A", 2, [20.0, 0.0, 0.0]),
        residue("B", 1, [0.0, 0.0, 3.0]),
        residue("B", 2, [20.0, 0.0, second_ligand_height]),
    ]
}

#[test]
fn a_model_identical_to_the_native_scores_one() {
    let native = two_sites(3.0);
    let result = score(&native, &native, options());
    assert!((result.fnat - 1.0).abs() < 1e-12);
    assert!(result.ligand_rmsd < 1e-5);
    assert!(result.interface_rmsd < 1e-5);
    assert!((result.score - 1.0).abs() < 1e-6, "score {}", result.score);
}

#[test]
fn fnat_counts_each_residue_pair_once() {
    // Many atom pairs touch within a residue pair, yet there are only two native
    // residue contacts. Lifting one ligand residue away loses exactly one of two.
    let result = score(&two_sites(33.0), &two_sites(3.0), options());
    assert!((result.fnat - 0.5).abs() < 1e-12, "fnat {}", result.fnat);
    // Four of the eight ligand backbone atoms moved 30 Å: sqrt(4 * 900 / 8).
    let expected = (4.0_f64 * 900.0 / 8.0).sqrt();
    assert!(
        (result.ligand_rmsd - expected).abs() < 1e-3,
        "ligand rmsd {}",
        result.ligand_rmsd
    );
}

#[test]
fn the_score_combines_the_components_as_published() {
    let result = score(&two_sites(33.0), &two_sites(3.0), options());
    let squash = |value: f64, scale: f64| 1.0 / (1.0 + (value / scale).powi(2));
    let expected =
        (result.fnat + squash(result.ligand_rmsd, 8.5) + squash(result.interface_rmsd, 1.5)) / 3.0;
    assert!((result.score - expected).abs() < 1e-12);
}

#[test]
fn a_rigidly_shifted_ligand_has_that_shift_as_its_ligand_rmsd() {
    let native = two_sites(3.0);
    let model = vec![
        residue("A", 1, [0.0, 0.0, 0.0]),
        residue("A", 2, [20.0, 0.0, 0.0]),
        residue("B", 1, [0.0, 0.0, 13.0]),
        residue("B", 2, [20.0, 0.0, 13.0]),
    ];
    let result = score(&model, &native, options());
    assert!(result.fnat.abs() < 1e-12);
    assert!(
        (result.ligand_rmsd - 10.0).abs() < 1e-4,
        "{}",
        result.ligand_rmsd
    );
}

#[test]
fn a_rigid_motion_of_the_whole_complex_does_not_change_the_score() {
    let native = two_sites(3.0);
    // Quarter turn about z then a translation: (x, y, z) -> (-y + 5, x - 7, z + 2).
    let mut text = String::new();
    let mut id = 0;
    for r in &native {
        for (element, name, o) in BACKBONE {
            id += 1;
            let (x, y, z) = (r.origin[0] + o[0], r.origin[1] + o[1], r.origin[2] + o[2]);
            let _ = writeln!(
                text,
                "ATOM {id} {element} {name} ALA {} {} {} {} {}",
                r.chain,
                r.number,
                -y + 5.0,
                x - 7.0,
                z + 2.0
            );
        }
    }
    let Ok(result) = dockq(
        &structure(&text),
        &structure(&body(&native)),
        "A",
        "B",
        options(),
    ) else {
        panic!("valid");
    };
    assert!((result.score - 1.0).abs() < 1e-5, "score {}", result.score);
}

#[test]
fn contact_and_interface_cutoffs_are_independent() {
    // Site one is a real contact; site two sits about 8 Å above its receptor
    // residue: an interface neighbour but not a contact.
    let native = two_sites(8.0);
    let mut model = two_sites(8.0);
    model[3] = residue("B", 2, [23.0, 0.0, 8.0]);
    let wide = score(&model, &native, options());
    let narrow = score(
        &model,
        &native,
        DockQOptions {
            interface_distance: 6.0,
            ..options()
        },
    );
    // Only the first site is a contact, and it is untouched.
    assert!((wide.fnat - 1.0).abs() < 1e-12);
    assert!((narrow.fnat - 1.0).abs() < 1e-12);
    assert!(wide.interface_rmsd > 0.1, "wide {}", wide.interface_rmsd);
    assert!(
        narrow.interface_rmsd < 1e-5,
        "narrow {}",
        narrow.interface_rmsd
    );
}

#[test]
fn hydrogens_do_not_make_contacts() {
    // The only atom close to the receptor is a ligand hydrogen.
    let mut native = vec![
        residue("A", 1, [0.0, 0.0, 0.0]),
        residue("B", 1, [0.0, 0.0, 9.0]),
    ];
    native[1].extra.push(("H", "H1", [0.0, 0.0, -8.0]));
    let outcome = dockq(
        &structure(&body(&native)),
        &structure(&body(&native)),
        "A",
        "B",
        options(),
    );
    assert!(matches!(outcome, Err(CompareError::NoComparablePairs)));
}

#[test]
fn cutoffs_must_be_positive_and_finite() {
    let native = two_sites(3.0);
    let model = structure(&body(&native));
    for bad in [0.0_f32, f32::NAN] {
        let outcome = dockq(
            &model,
            &model,
            "A",
            "B",
            DockQOptions {
                interface_distance: bad,
                ..options()
            },
        );
        assert!(matches!(outcome, Err(CompareError::InvalidDistanceCutoff)));
    }
}
