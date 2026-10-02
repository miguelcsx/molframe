use super::*;
use molframe_bench::{Sample, structure};

/// How many of the residues DSSP evaluated carry the same helix, strand or
/// coil class as the deposited annotation, and how many it evaluated.
fn agreement(sample: Sample) -> (usize, usize) {
    let structure = structure(sample);
    let assigned = assign_secondary_structure(&structure);
    let declared = structure.data().secondary_structure.to_vec();
    let mut compared = 0;
    let mut agreed = 0;
    for (declared, assigned) in declared.iter().zip(&assigned) {
        let class = |state: &SecondaryStructure| match state {
            SecondaryStructure::Helix => Some('H'),
            SecondaryStructure::Strand => Some('E'),
            SecondaryStructure::Coil | SecondaryStructure::Turn => Some('-'),
            SecondaryStructure::Unknown => None,
        };
        // A deposited annotation is silent about coil, so what it leaves
        // unlabelled on a residue DSSP evaluated is coil.
        let Some(assigned) = class(assigned) else {
            continue;
        };
        let declared = class(declared).unwrap_or('-');
        compared += 1;
        agreed += usize::from(declared == assigned);
    }
    (agreed, compared)
}

#[test]
fn non_polymer_residues_are_left_unknown() {
    let structure = structure(Sample::Small);
    let assigned = assign_secondary_structure(&structure);
    let unknown = assigned
        .iter()
        .filter(|state| **state == SecondaryStructure::Unknown)
        .count();
    // Ubiquitin's 76 residues are polymer; the remainder are waters.
    assert_eq!(assigned.len() - unknown, 76);
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
    let count = |wanted: SecondaryStructure| declared.iter().filter(|s| **s == wanted).count();
    assert_eq!(count(SecondaryStructure::Helix), 13 + 8);
    assert_eq!(count(SecondaryStructure::Strand), 4 + 4);
}
