use super::classify::classify;
use super::{Backbone, DsspOptions, hbond_energy, secondary_structure};
use molframe_chem::PolymerAtomRole;
use molframe_core::Presence;
use molframe_core::SecondaryStructure as Ss;
use molframe_core::annotation::{AnnotationColumn, AtomAnnotation};
use molframe_core::index::ResidueIndex;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;
use std::collections::BTreeSet;

fn options() -> DsspOptions {
    DsspOptions {
        electrostatic_prefactor: 332.0 * 0.42 * 0.20,
        hydrogen_bond_energy: -0.5,
        amide_hydrogen_distance: 1.0,
        minimum_sequence_separation: 2,
        helix_offset: 4,
        three_ten_offset: 3,
        pi_offset: 5,
        turn_offsets: 3..=5,
        bend_angle_degrees: 70.0,
    }
}

fn carbonyl(carbon: [f32; 3], oxygen: [f32; 3]) -> Backbone {
    Backbone {
        ca: None,
        nitrogen: None,
        carbon: Some(carbon),
        oxygen: Some(oxygen),
        hydrogen: None,
    }
}

fn amide(nitrogen: [f32; 3], hydrogen: [f32; 3]) -> Backbone {
    Backbone {
        ca: None,
        nitrogen: Some(nitrogen),
        carbon: None,
        oxygen: None,
        hydrogen: Some(hydrogen),
    }
}

#[test]
fn a_good_donor_acceptor_geometry_is_a_hydrogen_bond() {
    // O···N ≈ 2.9 Å with the hydrogen between them: a strong backbone bond.
    let energy = hbond_energy(
        &carbonyl([-1.23, 0.0, 0.0], [0.0, 0.0, 0.0]),
        &amide([2.9, 0.0, 0.0], [1.9, 0.0, 0.0]),
        options().electrostatic_prefactor,
    );
    assert!(energy < -0.5, "energy {energy}");
}

#[test]
fn a_distant_pair_is_not_a_hydrogen_bond() {
    let energy = hbond_energy(
        &carbonyl([-1.23, 0.0, 0.0], [0.0, 0.0, 0.0]),
        &amide([10.0, 0.0, 0.0], [9.0, 0.0, 0.0]),
        options().electrostatic_prefactor,
    );
    assert!(energy > -0.5, "energy {energy}");
}

/// Classifies residues without a Cα trace, so nothing can bend.
fn classified(bonds: &[(usize, usize)], evaluable: &[bool]) -> Vec<Ss> {
    let bonds: BTreeSet<(usize, usize)> = bonds.iter().copied().collect();
    classify(&bonds, evaluable, &vec![None; evaluable.len()], &options())
}

#[test]
fn consecutive_i_to_i_plus_four_bonds_make_a_helix() {
    let kinds = classified(&[(0, 4), (1, 5), (2, 6)], &[true; 8]);
    for kind in &kinds[1..=5] {
        assert_eq!(*kind, Ss::AlphaHelix);
    }
}

#[test]
fn a_lone_reciprocal_pair_is_an_isolated_bridge() {
    let kinds = classified(&[(2, 8), (8, 2)], &[true; 10]);
    assert_eq!(kinds[2], Ss::BetaBridge);
    assert_eq!(kinds[8], Ss::BetaBridge);
}

#[test]
fn a_lone_short_bond_makes_a_turn() {
    let kinds = classified(&[(0, 3)], &[true; 5]);
    assert_eq!(kinds[1], Ss::Turn);
    assert_eq!(kinds[2], Ss::Turn);
}

#[test]
fn a_short_peptide_runs_and_yields_one_record_per_residue() {
    let source = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 N N GLY A 1 0 0 0\n\
ATOM 2 C CA GLY A 1 1.5 0 0\n\
ATOM 3 C C GLY A 1 2 1 0\n\
ATOM 4 O O GLY A 1 3 1 0\n\
ATOM 5 N N GLY A 2 1.5 2 0\n\
ATOM 6 C CA GLY A 2 2 3 0\n\
ATOM 7 C C GLY A 2 3 3 0\n\
ATOM 8 O O GLY A 2 4 3 0\n";
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    let (structure, _): (Structure, _) = match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok(result) => result,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    };
    let structure = with_roles(&structure);
    let Ok(records) = secondary_structure(&structure, &options()) else {
        panic!("valid explicit DSSP definition");
    };
    assert_eq!(records.len(), 2);
    assert_eq!(records.kind(), [Ss::Unknown, Ss::Coil]);
    assert_eq!(
        records.placements().collect::<Vec<_>>(),
        vec![
            (ResidueIndex::new(0), Ss::Unknown),
            (ResidueIndex::new(1), Ss::Coil),
        ]
    );
}

#[test]
fn incomplete_backbones_are_unknown_instead_of_coil() {
    let kinds = classified(&[], &[true, false, true]);
    assert_eq!(kinds, [Ss::Coil, Ss::Unknown, Ss::Coil]);
}

#[test]
fn two_consecutive_three_turns_make_a_three_ten_helix() {
    let kinds = classified(&[(1, 4), (2, 5)], &[true; 8]);
    for kind in &kinds[2..=4] {
        assert_eq!(*kind, Ss::ThreeTenHelix);
    }
    assert_eq!(kinds[6], Ss::Coil);
}

#[test]
fn two_consecutive_five_turns_make_a_pi_helix() {
    let kinds = classified(&[(1, 6), (2, 7)], &[true; 10]);
    for kind in &kinds[2..=6] {
        assert_eq!(*kind, Ss::PiHelix);
    }
}

#[test]
fn bridges_in_register_form_a_ladder_of_strands() {
    let kinds = classified(&[(2, 9), (9, 2), (3, 8), (8, 3)], &[true; 12]);
    for residue in [2, 3, 8, 9] {
        assert_eq!(kinds[residue], Ss::Strand, "residue {residue}");
    }
}

#[test]
fn a_sharp_ca_angle_bends_and_the_threshold_is_an_option() {
    let trace = [
        Some([0.0, 0.0, 0.0]),
        Some([1.0, 0.0, 0.0]),
        Some([2.0, 0.0, 0.0]),
        Some([2.2, 1.0, 0.0]),
        Some([2.4, 2.0, 0.0]),
    ];
    let bonds = BTreeSet::new();
    let kinds = classify(&bonds, &[true; 5], &trace, &options());
    assert_eq!(kinds[2], Ss::Bend);
    assert_eq!(kinds[1], Ss::Coil);
    let lenient = DsspOptions {
        bend_angle_degrees: 85.0,
        ..options()
    };
    assert_eq!(classify(&bonds, &[true; 5], &trace, &lenient)[2], Ss::Coil);
}

fn with_roles(structure: &Structure) -> Structure {
    let entries: Vec<_> = structure
        .data()
        .atoms()
        .map(|atom| {
            let role = match atom.name() {
                Some("N") => PolymerAtomRole::PROTEIN_NITROGEN,
                Some("CA") => PolymerAtomRole::PROTEIN_ALPHA_CARBON,
                Some("C") => PolymerAtomRole::PROTEIN_CARBONYL_CARBON,
                Some("O") => PolymerAtomRole::PROTEIN_CARBONYL_OXYGEN,
                _ => PolymerAtomRole::UNKNOWN,
            };
            (role.code(), Presence::Present)
        })
        .collect();
    let mut data = structure.data().clone();
    let _ = data.annotations.insert(
        molframe_core::POLYMER_ATOM_ROLE_ANNOTATION,
        AtomAnnotation::Integer(AnnotationColumn::from_entries(entries).expect("small column")),
    );
    Structure::new(data)
}
