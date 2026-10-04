//! Pinned evidence: molframe against reference output of an independent library.
//!
//! Each test reads a JSON file under `tests/reference/` whose name carries the
//! library and its version, checks that the fixture it was generated from is the
//! fixture checked in here (by SHA-256), and compares the two within a tolerance
//! derived from the method, not from the observed difference.

use std::path::Path;

use molframe::ExecutionContext;
use molframe::chemistry::{RadiusSet, atom_radii};
use molframe_bench::{Sample, structure_from_cif};
use serde_json::Value;
use sha2::{Digest, Sha256};

fn reference(name: &str) -> Value {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/reference")
        .join(name);
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|error| panic!("{}: {error}", path.display()));
    serde_json::from_str(&text).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

fn numbers(value: &Value) -> Vec<f64> {
    let Some(items) = value.as_array() else {
        panic!("expected an array, found {value}");
    };
    items
        .iter()
        .map(|item| {
            item.as_f64()
                .unwrap_or_else(|| panic!("not a number: {item}"))
        })
        .collect()
}

fn indices(value: &Value) -> Vec<usize> {
    let Some(items) = value.as_array() else {
        panic!("expected an array, found {value}");
    };
    items
        .iter()
        .map(|item| {
            item.as_u64()
                .and_then(|index| usize::try_from(index).ok())
                .unwrap_or_else(|| panic!("not an index: {item}"))
        })
        .collect()
}

/// GW-013: Shrake–Rupley areas agree with biotite 1.7.1 atom by atom.
///
/// Both programs use molframe's Bondi radii, a 1.4 Å probe and 960 test points,
/// and neither lets water occlude. The point sets differ, so each atom's exposed
/// fraction is a different sample of the same surface; the standard error of a
/// fraction estimated from `n` points is at most `1/√n` of the sphere area
/// `4π(r + probe)²`. That bound is a property of the method, not a fitted number.
#[test]
fn gw_013_shrake_rupley_matches_biotite_per_atom() {
    let reference = reference("sasa_biotite_1.7.1.json");
    assert_eq!(reference["library"], "biotite");
    assert_eq!(reference["version"], "1.7.1");
    let probe = reference["probe"].as_f64().unwrap_or(f64::NAN);
    let points = reference["points"].as_u64().unwrap_or(0);
    let points = u16::try_from(points).unwrap_or(0);
    assert_eq!((probe, points), (1.4, 960));

    for (sample, key) in [(Sample::Tiny, "1crn"), (Sample::Small, "1ubq")] {
        let entry = &reference["entries"][key];
        let bytes = sample
            .cif()
            .unwrap_or_else(|| panic!("{key} ships no mmCIF"));
        assert_eq!(
            format!("{:x}", Sha256::digest(bytes)),
            entry["sha256"].as_str().unwrap_or_default(),
            "{key}: the fixture is not the one the reference was generated from"
        );

        let structure =
            structure_from_cif(sample).unwrap_or_else(|| panic!("{key} could not be read"));
        let kept = indices(&entry["atoms"]);
        let expected = numbers(&entry["areas"]);
        assert_eq!(kept.len(), expected.len());

        let all_positions = structure.positions();
        let all_radii = atom_radii(&structure, RadiusSet::Bondi);
        let positions: Vec<[f32; 3]> = kept.iter().map(|&atom| all_positions[atom]).collect();
        let radii: Vec<f32> = kept.iter().map(|&atom| all_radii[atom]).collect();
        let areas = molframe::surface::shrake_rupley(
            &positions,
            &radii,
            1.4,
            points,
            &ExecutionContext::default(),
        )
        .unwrap_or_else(|error| panic!("{key}: SASA failed: {error}"));
        assert_eq!(areas.len(), expected.len());

        let mut worst = 0.0_f64;
        for ((area, want), radius) in areas.iter().zip(&expected).zip(&radii) {
            let sphere = 4.0 * std::f64::consts::PI * (f64::from(*radius) + 1.4).powi(2);
            let tolerance = sphere / f64::from(points).sqrt();
            assert!(
                (area - want).abs() <= tolerance,
                "{key}: {area} against {want} (±{tolerance})"
            );
            worst = worst.max((area - want).abs() / sphere);
        }
        let total: f64 = areas.iter().sum();
        let want_total: f64 = expected.iter().sum();
        assert!(
            (total - want_total).abs() / want_total < 0.005,
            "{key}: total {total} against {want_total} (worst atom {worst:.4} of its sphere)"
        );
    }
}
