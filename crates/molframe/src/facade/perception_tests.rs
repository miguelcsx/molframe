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
HELIX    1   1 GLY A    1  GLY A    4  1                                   4\n\
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
fn file_secondary_structure_survives_automatic_assignment() {
    let structure = a_facade_read();
    assert_eq!(
        structure.engine().secondary_structure(),
        &[molframe_core::SecondaryStructure::AlphaHelix; 4]
    );
    assert_eq!(
        structure.secondary_source(),
        &[molframe_core::SecondarySource::File; 4]
    );
}

#[cfg(all(feature = "chemistry", feature = "spatial", feature = "mmcif"))]
#[test]
fn analysis_fills_only_the_residues_the_file_left_unassigned() {
    use molframe_core::{SecondarySource, SecondaryStructure};

    let Some(text) = molframe_bench::Sample::Tiny.cif() else {
        panic!("crambin ships as mmCIF")
    };
    let Ok((structure, _)) = read_bytes(text.to_vec(), Some("1crn.cif"), &ReadOptions::new())
    else {
        panic!("crambin reads")
    };
    let states = structure.secondary_structure();
    let sources = structure.secondary_source();
    assert_eq!(states.len(), sources.len());
    let from_file = sources
        .iter()
        .filter(|source| **source == SecondarySource::File)
        .count();
    // HELIX 7–19 and 23–30, SHEET 1–4 and 32–35.
    assert_eq!(from_file, 13 + 8 + 4 + 4);
    for (state, source) in states.iter().zip(sources) {
        match source {
            SecondarySource::File => {
                assert!(state.is_helix() || *state == SecondaryStructure::Strand);
            }
            SecondarySource::Dssp => assert_ne!(*state, SecondaryStructure::Unknown),
            SecondarySource::None => assert_eq!(*state, SecondaryStructure::Unknown),
            SecondarySource::CaOnly => panic!("crambin has full backbones"),
        }
    }
    assert!(sources.contains(&SecondarySource::Dssp));
}
