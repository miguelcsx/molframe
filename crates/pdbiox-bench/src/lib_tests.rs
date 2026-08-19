//! Confirms every embedded fixture parses through the shared loaders.

use super::*;

#[test]
fn every_sample_parses_from_binary_cif_with_atoms() {
    for sample in Sample::ALL {
        let structure = structure(sample);
        assert!(
            structure.atom_count() > 0,
            "{} parsed to zero atoms",
            sample.label()
        );
    }
}

#[test]
fn legacy_pdb_and_mmcif_agree_on_atom_count() {
    for sample in Sample::MODELS {
        let reference = structure(sample).atom_count();
        if let Some(from_pdb) = structure_from_pdb(sample) {
            assert_eq!(from_pdb.atom_count(), reference, "{}", sample.label());
        }
        if let Some(from_cif) = structure_from_cif(sample) {
            assert_eq!(from_cif.atom_count(), reference, "{}", sample.label());
        }
    }
}

#[test]
fn the_ensemble_has_several_models() {
    let structure = structure(Sample::Ensemble);
    assert!(
        structure.model_count() > 1,
        "ensemble collapsed to one model"
    );
    assert!(model_coordinates(&structure, 1).is_some());
    assert!(
        structure.atom_count() > 100,
        "ensemble fixture parsed with too few atoms"
    );
}

#[test]
fn the_large_sample_decompresses_from_gzip() {
    let buffer = input_gzip(large_cif_gz());
    let parsed = pdbiox_cif::read(&buffer, &pdbiox_core::io::ReadOptions::new());
    match parsed {
        Ok((structure, _)) => assert!(structure.atom_count() > 10_000),
        Err(findings) => panic!("gzip mmCIF failed: {findings:?}"),
    }
}

#[test]
fn coordinate_helpers_track_the_structure() {
    let structure = structure(Sample::Small);
    let coordinates = coordinates(&structure);
    let Ok(coordinate_count) = u32::try_from(coordinates.len()) else {
        panic!("coordinate fixture is too large for the structure index")
    };
    assert_eq!(coordinate_count, structure.atom_count());

    let moved = perturbed(&coordinates, 0.5);
    assert_eq!(moved.len(), coordinates.len());
    for (original, shifted) in coordinates.iter().zip(&moved) {
        for axis in 0..3 {
            assert!((original[axis] - shifted[axis]).abs() <= 0.5);
        }
    }

    assert_eq!(
        uniform_radii(coordinates.len(), 1.7).len(),
        coordinates.len()
    );
}

#[test]
fn component_dictionary_fixtures_are_present() {
    assert!(!ccd_hem().is_empty());
    assert!(!ccd_atp().is_empty());
}
