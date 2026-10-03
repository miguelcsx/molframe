use super::*;
use molframe_bench::{Sample, structure};

/// The three-state class DSSP agreement is reported in: H (any helix), E
/// (strand or bridge) or coil; `None` for a residue nothing assigned.
fn three_state(state: SecondaryStructure) -> Option<char> {
    match state {
        SecondaryStructure::Unknown => None,
        state if state.is_helix() => Some('H'),
        SecondaryStructure::Strand | SecondaryStructure::BetaBridge => Some('E'),
        _ => Some('-'),
    }
}

/// How many of the residues DSSP evaluated carry the same three-state class
/// as the deposited annotation, and how many it evaluated.
fn agreement(sample: Sample) -> (usize, usize) {
    let structure = structure(sample);
    let assigned = assign_secondary_structure(&structure);
    let declared = structure.data().secondary_structure.to_vec();
    let mut compared = 0;
    let mut agreed = 0;
    for (declared, assigned) in declared.iter().zip(&assigned) {
        // A deposited annotation is silent about coil, so what it leaves
        // unlabelled on a residue DSSP evaluated is coil.
        let Some(assigned) = three_state(assigned.state) else {
            continue;
        };
        let declared = three_state(*declared).unwrap_or('-');
        compared += 1;
        agreed += usize::from(declared == assigned);
    }
    (agreed, compared)
}

#[test]
fn non_polymer_residues_are_left_unknown_with_no_source() {
    let structure = structure(Sample::Small);
    let assigned = assign_secondary_structure(&structure);
    let unknown = assigned
        .iter()
        .filter(|assignment| assignment.state == SecondaryStructure::Unknown)
        .count();
    // Ubiquitin's 76 residues are polymer; the remainder are waters.
    assert_eq!(assigned.len() - unknown, 76);
    assert!(assigned.iter().all(|assignment| {
        (assignment.state == SecondaryStructure::Unknown)
            == (assignment.source == SecondarySource::None)
    }));
    assert!(
        assigned
            .iter()
            .filter(|assignment| assignment.state != SecondaryStructure::Unknown)
            .all(|assignment| assignment.source == SecondarySource::Dssp)
    );
}

#[test]
fn dssp_agrees_with_the_deposited_annotation_on_most_residues() {
    for (sample, minimum) in [
        (Sample::Tiny, 0.75),
        (Sample::Small, 0.80),
        (Sample::Medium, 0.80),
        (Sample::Large, 0.95),
    ] {
        let (agreed, compared) = agreement(sample);
        assert!(compared > 0);
        let fraction = f64::from(u32::try_from(agreed).unwrap_or(u32::MAX))
            / f64::from(u32::try_from(compared).unwrap_or(u32::MAX));
        assert!(
            fraction >= minimum,
            "{agreed} of {compared} residues agree, below {minimum}"
        );
    }
}

/// Crambin's deposited HELIX 7–19 and 23–30 and SHEET 1–4 and 32–35 are the
/// ranges the binary reader must hand over, not leave for DSSP to rediscover.
#[test]
fn the_binary_reader_keeps_the_deposited_helix_and_sheet_ranges() {
    let declared = structure(Sample::Tiny).data().secondary_structure.to_vec();
    assert_eq!(declared.iter().filter(|s| s.is_helix()).count(), 13 + 8);
    assert_eq!(
        declared
            .iter()
            .filter(|s| **s == SecondaryStructure::Strand)
            .count(),
        4 + 4
    );
}

/// Classifies `residues` coil residues from hand-written hydrogen bonds.
fn classified(
    residues: usize,
    bonds: &[(usize, usize)],
    pairs: &[(usize, usize)],
) -> Vec<SecondaryStructure> {
    let backbones = vec![Backbone::default(); residues];
    let bonds: Bonds = bonds.iter().copied().collect();
    let mut states = vec![SecondaryStructure::Coil; residues];
    classify(&backbones, &bonds, pairs, &mut states);
    states
}

#[test]
fn two_consecutive_three_turns_make_a_three_ten_helix() {
    let states = classified(8, &[(1, 4), (2, 5)], &[]);
    assert_eq!(states[2..=4], [SecondaryStructure::ThreeTenHelix; 3]);
    assert_eq!(states[1], SecondaryStructure::Coil);
    assert_eq!(states[5], SecondaryStructure::Coil);
}

#[test]
fn two_consecutive_five_turns_make_a_pi_helix() {
    let states = classified(10, &[(1, 6), (2, 7)], &[]);
    assert!(
        states[2..=6]
            .iter()
            .all(|state| *state == SecondaryStructure::PiHelix)
    );
}

#[test]
fn a_single_bridge_outside_a_ladder_is_a_beta_bridge() {
    let states = classified(12, &[(2, 9), (9, 2)], &[(2, 9)]);
    assert_eq!(states[2], SecondaryStructure::BetaBridge);
    assert_eq!(states[9], SecondaryStructure::BetaBridge);
    assert_eq!(states[3], SecondaryStructure::Coil);
}

#[test]
fn a_sharp_ca_angle_is_a_bend_and_a_straight_one_is_not() {
    let trace = |positions: [[f32; 3]; 5]| {
        positions
            .map(|position| Backbone {
                ca: Some(position),
                carbon: Some(position),
                nitrogen: Some(position),
                ..Backbone::default()
            })
            .to_vec()
    };
    let bent = trace([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [2.2, 1.0, 0.0],
        [2.4, 2.0, 0.0],
    ]);
    assert!(bends(&bent, 2));
    assert!(!bends(&bent, 1), "residue 1 has no Cα two before it");
    let straight = trace([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [4.0, 0.0, 0.0],
    ]);
    assert!(!bends(&straight, 2));

    let mut states = vec![SecondaryStructure::Coil; 5];
    classify(&bent, &Bonds::default(), &[], &mut states);
    assert_eq!(states[2], SecondaryStructure::Bend);
}
