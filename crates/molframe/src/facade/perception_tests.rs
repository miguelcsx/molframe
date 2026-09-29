//! Default perception on a facade read: bonds and secondary structure.

use super::*;

#[cfg(feature = "pdb")]
fn a_facade_read() -> Structure {
    let Ok((structure, _)) = read_bytes(
        DEFAULT_PERCEPTION.as_bytes().to_vec(),
        Some("perception.pdb"),
        &ReadOptions::new(),
    ) else {
        panic!("perception fixture must parse")
    };
    structure
}

/// One glycine whose backbone bonds and a HELIX record must survive the
/// automatic fallback assignment of a default facade read.
const DEFAULT_PERCEPTION: &str = "\
HELIX    1   1 GLY A   1  GLY A   4  1                                  4\n\
ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00 10.00           N\n\
ATOM      2  CA  GLY A   1       1.450   0.000   0.000  1.00 10.00           C\n\
ATOM      3  C   GLY A   1       2.900   0.000   0.000  1.00 10.00           C\n\
ATOM      4  O   GLY A   1       4.100   0.000   0.000  1.00 10.00           O\n\
ATOM      5  CA  GLY A   2       2.100   2.000   0.000  1.00 10.00           C\n\
ATOM      6  CA  GLY A   3       4.200   2.000   0.000  1.00 10.00           C\n\
ATOM      7  CA  GLY A   4       6.300   2.000   0.000  1.00 10.00           C\n\
END\n\
";

#[cfg(all(feature = "chemistry", feature = "spatial", feature = "pdb"))]
#[test]
fn default_read_adds_distance_bonds_with_inferred_provenance() {
    let bonds: Vec<_> = a_facade_read().engine().data().bonds.iter().collect();
    assert!(bonds.len() >= 3, "bonds: {bonds:?}");
    assert!(
        bonds
            .iter()
            .all(|bond| bond.provenance == molframe_core::BondProvenance::InferredDistance)
    );
}

#[cfg(all(feature = "chemistry", feature = "spatial", feature = "pdb"))]
#[test]
fn default_read_fills_only_the_secondary_structure_the_file_left_unknown() {
    assert_eq!(
        a_facade_read().engine().secondary_structure(),
        &[molframe_core::SecondaryStructure::Helix; 4]
    );
}

#[cfg(all(feature = "chemistry", feature = "spatial", feature = "pdb"))]
#[test]
fn file_secondary_structure_survives_automatic_fallback_assignment() {
    assert_eq!(
        a_facade_read().engine().secondary_structure(),
        &[molframe_core::SecondaryStructure::Helix; 4]
    );
}
