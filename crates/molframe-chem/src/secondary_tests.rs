use super::classify::classify;
use super::geometry::{amide_hydrogen, bends, bond_energy};
use super::*;
use DsspBackbone as Backbone;
use molframe_bench::{Sample, structure};
use num_traits::ToPrimitive;

fn continuous_backbones(count: usize) -> Vec<Backbone> {
    (0..count)
        .map(|i| {
            let x = i.to_f32().expect("small test index");
            Backbone {
                ca: Some([x, 0.0, 0.0]),
                carbon: Some([x, 0.0, 0.0]),
                nitrogen: Some([x, 0.0, 0.0]),
                oxygen: Some([x, 1.0, 0.0]),
                ..Backbone::default()
            }
        })
        .collect()
}
fn classified(count: usize, bonds: &[(usize, usize)], pairs: &[(usize, usize)]) -> Vec<Ss> {
    let mut states = vec![Ss::Coil; count];
    classify(
        &continuous_backbones(count),
        &bonds.iter().copied().collect(),
        pairs,
        &mut states,
        &DsspOptions::default(),
    );
    states
}

#[test]
fn non_polymer_residues_are_left_unknown_with_no_source() {
    let assigned = assign_secondary_structure(&structure(Sample::Small));
    assert_eq!(
        assigned.iter().filter(|a| a.state != Ss::Unknown).count(),
        76
    );
    assert!(
        assigned
            .iter()
            .all(|a| (a.state == Ss::Unknown) == (a.source == SecondarySource::None))
    );
}
#[test]
fn the_binary_reader_keeps_deposited_helix_and_sheet_ranges() {
    let declared = structure(Sample::Tiny).data().secondary_structure.to_vec();
    assert_eq!(declared.iter().filter(|s| s.is_helix()).count(), 21);
    assert_eq!(declared.iter().filter(|s| **s == Ss::Strand).count(), 8);
}
#[test]
fn consecutive_turns_make_alpha_three_ten_and_pi_helices() {
    for (size, kind) in [
        (3, Ss::ThreeTenHelix),
        (4, Ss::AlphaHelix),
        (5, Ss::PiHelix),
    ] {
        let states = classified(10, &[(1, 1 + size), (2, 2 + size)], &[]);
        assert!(states[2..2 + size].iter().all(|s| *s == kind));
    }
    assert_eq!(
        classified(8, &[(1, 5)], &[])[2],
        Ss::Turn,
        "one turn is not a helix"
    );
}
#[test]
fn an_occupied_position_rejects_the_entire_three_ten_stretch() {
    let states = classified(12, &[(1, 4), (2, 5), (4, 9), (9, 4)], &[(4, 9)]);
    assert_eq!(states[2], Ss::Turn);
    assert_eq!(states[3], Ss::Turn);
    assert_eq!(states[4], Ss::BetaBridge);
}
#[test]
fn pi_helices_can_replace_alpha_but_not_a_sheet() {
    let states = classified(12, &[(1, 5), (2, 6), (1, 6), (2, 7)], &[]);
    assert!(states[2..7].iter().all(|s| *s == Ss::PiHelix));
}
#[test]
fn cross_chain_ladders_and_beta_bulges_fill_their_whole_spans() {
    let mut backbones = continuous_backbones(20);
    for b in &mut backbones[10..] {
        b.chain = 1;
    }
    let bonds = [(2, 17), (17, 2), (3, 16), (16, 3), (5, 14), (14, 5)];
    let mut states = vec![Ss::Coil; 20];
    classify(
        &backbones,
        &bonds.into_iter().collect(),
        &[(2, 17), (3, 16), (5, 14)],
        &mut states,
        &DsspOptions::default(),
    );
    assert!(states[2..=5].iter().all(|s| *s == Ss::Strand));
    assert!(states[14..=17].iter().all(|s| *s == Ss::Strand));
    assert_eq!(states[4], Ss::Strand, "bulge interior has no direct bridge");
}
#[test]
fn peptide_breaks_and_chain_boundaries_block_local_patterns() {
    let mut backbones = continuous_backbones(10);
    backbones[3].chain = 1;
    let mut states = vec![Ss::Coil; 10];
    classify(
        &backbones,
        &[(1, 5), (2, 6)].into_iter().collect(),
        &[],
        &mut states,
        &DsspOptions::default(),
    );
    assert!(states.iter().all(|s| *s == Ss::Coil));
    backbones[3].chain = 0;
    backbones[3].nitrogen = Some([100.0, 0.0, 0.0]);
    classify(
        &backbones,
        &[(1, 5), (2, 6)].into_iter().collect(),
        &[],
        &mut states,
        &DsspOptions::default(),
    );
    assert!(states.iter().all(|s| *s == Ss::Coil));
}
#[test]
fn proline_cannot_donate_and_oxygen_is_required_for_evaluation() {
    let mut backbones = continuous_backbones(3);
    backbones[1].proline = true;
    assert_eq!(amide_hydrogen(&backbones, 1, 1.0), None);
    assert!(bond_energy(&backbones, 0, 1, &DsspOptions::default()).abs() < f64::EPSILON);
    backbones[2].oxygen = None;
    assert!(!backbones[2].is_evaluable());
    assert_eq!(
        dssp_from_backbones(&backbones, &DsspOptions::default()).expect("valid options")[2],
        Ss::Unknown
    );
}

#[test]
fn native_options_reject_invalid_numbers_and_ca_fallback_can_name_strands() {
    let invalid = DsspOptions {
        hydrogen_bond_energy: f64::NAN,
        ..DsspOptions::default()
    };
    assert_eq!(dssp_from_backbones(&[], &invalid), Err(InvalidDsspOptions));
    let mut backbones = continuous_backbones(4);
    for (i, b) in backbones.iter_mut().enumerate() {
        b.ca = Some([i.to_f32().expect("small index") * 3.5, 0.0, 0.0]);
        b.carbon = None;
        b.nitrogen = None;
        b.oxygen = None;
    }
    let mut states = vec![Ss::Coil; 4];
    assign_ca_trace(&backbones, &mut states);
    assert_eq!(states, [Ss::Strand; 4]);
    backbones[2].ca = Some([100.0, 0.0, 0.0]);
    states.fill(Ss::Coil);
    assign_ca_trace(&backbones, &mut states);
    assert_eq!(states, [Ss::Coil; 4]);
}
#[test]
fn a_sharp_ca_direction_is_a_bend() {
    let positions = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [2.2, 1.0, 0.0],
        [2.4, 2.0, 0.0],
    ];
    let backbones = positions.map(|p| Backbone {
        ca: Some(p),
        carbon: Some(p),
        nitrogen: Some(p),
        oxygen: Some(p),
        ..Backbone::default()
    });
    assert!(bends(&backbones, 2, 70.0));
    assert!(!bends(&backbones, 2, 85.0));
}

// 7QPD author B residues 2–6, RCSB CC0 coordinates. mkdssp 4.5 assigns P
// to residues 3–5, including ARG and THR: the rule is not proline-specific.
fn pp_backbones() -> Vec<Backbone> {
    [
        (
            [194.468, 175.337, 185.460],
            [195.493, 174.908, 184.509],
            [195.445, 175.804, 183.271],
            [194.784, 175.518, 182.271],
        ),
        (
            [196.169, 176.915, 183.356],
            [196.234, 177.844, 182.237],
            [197.106, 177.276, 181.124],
            [198.131, 176.641, 181.383],
        ),
        (
            [196.689, 177.501, 179.882],
            [197.414, 176.984, 178.735],
            [198.644, 177.841, 178.446],
            [198.632, 179.053, 178.672],
        ),
        (
            [199.719, 177.235, 177.946],
            [200.875, 178.032, 177.526],
            [200.549, 178.856, 176.294],
            [199.734, 178.460, 175.458],
        ),
        (
            [201.195, 180.011, 176.185],
            [201.002, 180.916, 175.057],
            [202.335, 181.057, 174.333],
            [203.148, 181.921, 174.671],
        ),
    ]
    .map(|(n, ca, c, o)| Backbone {
        nitrogen: Some(n),
        ca: Some(ca),
        carbon: Some(c),
        oxygen: Some(o),
        ..Backbone::default()
    })
    .to_vec()
}
#[test]
fn pp_requires_three_consecutive_phi_psi_windows_and_preserves_turns_and_bends() {
    let b = pp_backbones();
    let mut states = vec![Ss::Coil; 5];
    classify(
        &b,
        &Bonds::default(),
        &[],
        &mut states,
        &DsspOptions::default(),
    );
    assert_eq!(states[1..4], [Ss::PolyProline; 3]);
    let mut states = vec![Ss::Coil; 5];
    states[1] = Ss::Turn;
    states[2] = Ss::Bend;
    classify(
        &b,
        &Bonds::default(),
        &[],
        &mut states,
        &DsspOptions::default(),
    );
    assert_eq!(states[1..4], [Ss::Turn, Ss::Bend, Ss::PolyProline]);
    let mut states = vec![Ss::Coil; 4];
    classify(
        &b[..4],
        &Bonds::default(),
        &[],
        &mut states,
        &DsspOptions::default(),
    );
    assert!(!states.contains(&Ss::PolyProline));
    let mut broken = b;
    broken[2].chain = 1;
    let mut states = vec![Ss::Coil; 5];
    classify(
        &broken,
        &Bonds::default(),
        &[],
        &mut states,
        &DsspOptions::default(),
    );
    assert!(!states.contains(&Ss::PolyProline));
}

#[test]
fn every_1aon_exact_state_and_author_identity_matches_mkdssp_4_5() {
    let source = structure(Sample::Large);
    let assignments = assign_secondary_structure(&source);
    let mut rows = Vec::new();
    for chain in source.data().chains() {
        for residue in chain.residues() {
            let kind = assignments[residue.index().as_usize()].state;
            let code = match kind {
                Ss::Unknown => continue,
                Ss::AlphaHelix => 'H',
                Ss::ThreeTenHelix => 'G',
                Ss::PiHelix => 'I',
                Ss::PolyProline => 'P',
                Ss::Strand => 'E',
                Ss::BetaBridge => 'B',
                Ss::Turn => 'T',
                Ss::Bend => 'S',
                Ss::Coil => 'C',
                Ss::OtherHelix => panic!("not a native state"),
            };
            rows.push((
                chain.auth_label().expect("author chain"),
                residue.auth_seq_id().expect("author residue"),
                residue.ins_code().unwrap_or(""),
                code,
            ));
        }
    }
    rows.sort_unstable();
    assert_eq!(rows.len(), 8015);
    let mut digest = 14_695_981_039_346_656_037_u64;
    for (chain, seq, ins, code) in rows {
        for byte in format!("{chain}|{seq}|{ins}|{code}\n").bytes() {
            digest = (digest ^ u64::from(byte)).wrapping_mul(1_099_511_628_211);
        }
    }
    assert_eq!(
        digest, 0xb127_15e8_a8e2_14fd,
        "mkdssp 4.5.0 exact-state and identity digest"
    );
}
