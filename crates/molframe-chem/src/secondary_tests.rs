use super::*;
use molframe_bench::{Sample, structure};

fn counts(sample: Sample) -> [usize; 5] {
    let states = assign_secondary_structure(&structure(sample));
    let count =
        |wanted: SecondaryStructure| states.iter().filter(|state| **state == wanted).count();
    [
        states.len(),
        count(SecondaryStructure::Helix),
        count(SecondaryStructure::Strand),
        count(SecondaryStructure::Turn),
        count(SecondaryStructure::Coil),
    ]
}

/// Residues, helix, strand, turn and coil, as the grid-bounded pass assigned
/// them before the neighbour search moved onto the shared cell grid.
#[test]
fn assignment_counts_are_unchanged_on_the_embedded_structures() {
    assert_eq!(counts(Sample::Tiny), [46, 28, 4, 0, 14]);
    assert_eq!(counts(Sample::Small), [134, 29, 14, 10, 23]);
    assert_eq!(counts(Sample::Medium), [801, 489, 0, 32, 53]);
    assert_eq!(counts(Sample::Large), [8029, 4591, 586, 423, 2415]);
}
