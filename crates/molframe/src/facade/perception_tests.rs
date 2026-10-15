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

#[cfg(all(feature = "chemistry", feature = "spatial", feature = "pdb"))]
#[test]
fn skipping_perception_returns_the_structure_as_the_file_states_it() {
    let Ok((structure, _)) = read_bytes(
        DEFAULT_PERCEPTION.as_bytes().to_vec(),
        Some("perception.pdb"),
        &ReadOptions::new().perceive(false),
    ) else {
        panic!("perception fixture must parse")
    };
    assert_eq!(structure.engine().data().bonds.iter().count(), 0);
    assert_eq!(
        structure.engine().secondary_structure(),
        &[molframe_core::SecondaryStructure::AlphaHelix; 4]
    );
}

/// The perception fixture read without perception, then enriched against
/// `store` alone: no environment, no global store.
#[cfg(all(feature = "chemistry", feature = "spatial", feature = "pdb"))]
fn enriched_with(
    store: &molframe_chem::TableStore,
) -> (
    Vec<(molframe_core::AtomIndex, molframe_core::AtomIndex)>,
    Vec<Diagnostic>,
) {
    let mut options = ReadOptions::new();
    options.perception = molframe_core::io::Perception::Skip;
    let Ok((raw, _)) = read_bytes(
        DEFAULT_PERCEPTION.as_bytes().to_vec(),
        Some("perception.pdb"),
        &options,
    ) else {
        panic!("perception fixture must parse")
    };
    let Ok((structure, findings)) = super::enrich::enrich_read_with(
        Ok((raw.engine().clone(), Vec::new())),
        &ReadOptions::new(),
        &ExecutionContext::default(),
        store,
    ) else {
        panic!("the fixture enriches")
    };
    let bonds = structure
        .data()
        .bonds
        .iter()
        .map(|bond| (bond.atom_a, bond.atom_b))
        .collect();
    (bonds, findings)
}

/// A store over an empty directory with no profile file, so only what the
/// store holds in memory or compiles in can answer.
#[cfg(all(feature = "chemistry", feature = "spatial", feature = "pdb"))]
fn empty_store(name: &str) -> molframe_chem::TableStore {
    let dir = std::env::temp_dir().join(format!("molframe-facade-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    let _ = std::fs::create_dir_all(&dir);
    molframe_chem::TableStore::new(dir.clone()).with_profile(dir.join("profile.txt"))
}

#[cfg(all(feature = "chemistry", feature = "spatial", feature = "pdb"))]
fn shipped_reach() -> Vec<u8> {
    let path = format!(
        "{}/../../data/macromolecular-reach.reach",
        env!("CARGO_MANIFEST_DIR")
    );
    match std::fs::read(&path) {
        Ok(bytes) => bytes,
        Err(error) => panic!("{path}: {error}"),
    }
}

#[cfg(all(feature = "chemistry", feature = "spatial", feature = "pdb"))]
#[test]
fn a_reach_table_registered_in_memory_perceives_bonds_with_no_files() {
    let store = empty_store("memory");
    store.insert_profile("bond-reach: held\n");
    store.insert_table::<molframe_chem::BondReach>("held", shipped_reach());
    let (bonds, findings) = enriched_with(&store);
    assert!(bonds.len() >= 3, "bonds: {bonds:?}");
    assert!(findings.iter().all(|finding| finding.code() != Code::W3207));
}

#[cfg(all(
    feature = "chemistry",
    feature = "spatial",
    feature = "pdb",
    not(feature = "embedded-tables")
))]
#[test]
fn a_store_with_nothing_to_offer_skips_bonds_and_says_so() {
    let (bonds, findings) = enriched_with(&empty_store("nothing"));
    assert!(bonds.is_empty());
    assert!(findings.iter().any(|finding| finding.code() == Code::W3207));
}

#[cfg(all(feature = "embedded-tables", feature = "spatial", feature = "pdb"))]
#[test]
fn embedded_tables_perceive_bonds_with_no_profile_or_cache() {
    let store = empty_store("embedded");
    assert!(store.has_embedded_tables());
    let (bonds, findings) = enriched_with(&store);
    assert!(bonds.len() >= 3, "bonds: {bonds:?}");
    assert!(findings.iter().all(|finding| finding.code() != Code::W3207));
}

#[cfg(all(feature = "embedded-tables", feature = "spatial", feature = "pdb"))]
#[test]
fn embedded_tables_perceive_the_bonds_the_cache_file_does() {
    let embedded = enriched_with(&empty_store("embedded-same")).0;
    let store = empty_store("cache-same");
    let _ = std::fs::write(store.directory().join("profile.txt"), "bond-reach: lab\n");
    let _ = std::fs::write(store.directory().join("lab.reach"), shipped_reach());
    assert_eq!(embedded, enriched_with(&store).0);
}
