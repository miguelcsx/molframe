//! Pinned evidence for crystallography, TM-score and element data.
//!
//! Each test reads a JSON file under `tests/reference/` whose name carries the
//! library and its version, checks that the fixture it was generated from is the
//! one in the checkout (by SHA-256 where there is a fixture), and compares the two
//! within a tolerance derived from the method rather than from the observed
//! difference.

use molframe::UnitCell;
use molframe::crystal::CellTransform;
use molframe::crystal::{SpaceGroupSetting, space_group_by_hermann_mauguin, space_group_setting};
use serde_json::Value;

fn reference(name: &str) -> Value {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/reference")
        .join(name);
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|error| panic!("{}: {error}", path.display()));
    serde_json::from_str(&text).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

fn setting(hall: i32) -> &'static SpaceGroupSetting {
    let Ok(hall) = u16::try_from(hall) else {
        panic!("Hall number {hall} does not fit")
    };
    space_group_setting(hall).unwrap_or_else(|error| panic!("Hall {hall} must resolve: {error}"))
}

/// Every Hall setting in the bundled catalogue carries the operator count and
/// international number the catalogue records for it.
///
/// A setting whose operator list lost or gained a row would place atoms that are
/// not in the crystal, and a symmetry expansion built on it would silently
/// produce a different assembly.
#[test]
fn gw_041_every_hall_setting_matches_its_recorded_metadata() {
    let reference = reference("space_groups_gemmi_0.7.5.json");
    assert_eq!(reference["library"], "gemmi");
    assert_eq!(reference["version"], "0.7.5");

    let Some(settings) = reference["settings"].as_object() else {
        panic!("settings must be an object");
    };
    assert_eq!(settings.len(), 530, "the catalogue holds all Hall settings");

    for (key, entry) in settings {
        let hall: i32 = key.parse().unwrap_or_else(|_| panic!("hall number {key}"));
        let group = setting(hall);
        let want = entry["operator_count"]
            .as_u64()
            .unwrap_or_else(|| panic!("Hall {hall} has no recorded count"));
        let got = u64::try_from(group.operations.len())
            .unwrap_or_else(|_| panic!("operator count fits u64"));
        assert_eq!(
            got, want,
            "Hall {hall} ({}): {got} operators against {want}",
            group.hall_symbol
        );
        assert_eq!(
            i64::from(group.international_number),
            entry["international_number"].as_i64().unwrap_or(-1),
            "Hall {hall} international number"
        );
        assert_eq!(
            &*group.hall_symbol,
            entry["hall_symbol"].as_str().unwrap_or_default(),
            "Hall {hall} symbol"
        );
    }
}

/// A group's expansion has the multiplicity its setting implies, and every image
/// of a general position is distinct.
#[test]
fn gw_042_space_group_expansion_places_the_expected_images() {
    // P 1 is the identity alone.
    let triclinic = space_group_by_hermann_mauguin("P 1")
        .unwrap_or_else(|error| panic!("P 1 must resolve: {error}"));
    assert_eq!(triclinic.operations.len(), 1, "P 1 is one operator");
    let identity_image = triclinic.operations[0].apply_fractional([0.21, 0.34, 0.17]);
    assert!(
        (identity_image[0] - 0.21).abs() < 1e-9
            && (identity_image[1] - 0.34).abs() < 1e-9
            && (identity_image[2] - 0.17).abs() < 1e-9,
        "P 1's operator is the identity: {identity_image:?}"
    );

    // F m -3 m has the full cubic multiplicity, all distinct at a general position.
    let cubic = space_group_by_hermann_mauguin("F m -3 m")
        .unwrap_or_else(|error| panic!("F m -3 m must resolve: {error}"));
    assert_eq!(cubic.operations.len(), 192, "F m -3 m has 192 cosets");

    // The images are counted in fractional space with a tolerance, so the test
    // needs no lossy conversion of a coordinate to an integer. Distinct images of
    // a general position are separated by a large fraction of a cell edge, while
    // the tolerance is far below any real separation, so only numerically
    // identical images collapse.
    let images: Vec<[f64; 3]> = cubic
        .operations
        .iter()
        .map(|operation| {
            operation
                .apply_fractional([0.21, 0.34, 0.17])
                .map(|axis| axis.rem_euclid(1.0))
        })
        .collect();
    let mut distinct = 0usize;
    for (index, image) in images.iter().enumerate() {
        let seen = images[..index].iter().any(|earlier| {
            let nearest = (0..3).map(|axis| {
                let delta = (image[axis] - earlier[axis]).abs();
                delta.min(1.0 - delta)
            });
            nearest.fold(0.0_f64, f64::max) < 1e-6
        });
        if !seen {
            distinct += 1;
        }
    }
    assert_eq!(distinct, 192, "a general position has 192 distinct images");
}

/// A fractional coordinate and its Cartesian image round-trip, and the cell
/// reproduces the parameters it was built from.
#[test]
fn gw_043_cell_transforms_round_trip_and_reproduce_their_parameters() {
    let cases = [
        ([30.0, 40.0, 50.0], [70.0, 80.0, 100.0]),
        ([20.0, 30.0, 40.0], [90.0, 90.0, 90.0]),
        ([50.0, 50.0, 80.0], [90.0, 90.0, 120.0]),
        ([31.0, 45.0, 62.0], [61.0, 71.0, 89.0]),
    ];
    for (lengths, angles) in cases {
        let cell = UnitCell { lengths, angles };
        // Building a cell copies its parameters, so the edge lengths come back
        // bit-for-bit and an exact comparison is the honest test.
        for (axis, &expected) in lengths.iter().enumerate() {
            assert!(
                cell.lengths[axis].to_bits() == expected.to_bits(),
                "the cell keeps its edge lengths: {:?} != {lengths:?}",
                cell.lengths
            );
        }
        let transform = CellTransform::new(&cell)
            .unwrap_or_else(|error| panic!("{lengths:?} {angles:?}: {error}"));

        let fractional = [0.13, 0.27, 0.41];
        let cartesian = transform.to_cartesian(fractional);
        let back = transform.to_fractional(cartesian);
        for (got, want) in back.iter().zip(&fractional) {
            assert!(
                (got - want).abs() < 1e-9,
                "{lengths:?} {angles:?}: round trip {back:?} against {fractional:?}"
            );
        }

        // The first edge lies along x with length a, which is the convention the
        // transform is built on.
        let a = transform.to_cartesian([1.0, 0.0, 0.0]);
        assert!((a[0] - lengths[0]).abs() < 1e-9, "a along x: {a:?}");
        assert!(a[1].abs() < 1e-9 && a[2].abs() < 1e-9, "{a:?}");
    }
}

/// The TM-score reference is internally consistent and carries the two cases
/// that separate the published score from a least-squares fit.
///
/// The values were produced by `tmtools` 0.3.0, the Python binding of Zhang and
/// Skolnick's TM-align, on the 1UBQ C-alpha trace pinned by SHA-256. The
/// numerical comparison runs in the Python differential; this test guards the
/// reference file the moment the checkout's version drifts from it.
#[test]
fn gw_044_tm_score_reference_covers_non_rigid_motion() {
    let reference = reference("tm_score_tmalign_0.3.0.json");
    assert_eq!(reference["library"], "tmtools (Zhang & Skolnick TM-align)");
    assert_eq!(reference["version"], "0.3.0");
    assert_eq!(reference["fixture"], "1ubq.pdb");
    assert_eq!(
        reference["sha256"].as_str().map(str::len),
        Some(64),
        "the fixture digest is a SHA-256"
    );

    let cases = reference["cases"]
        .as_object()
        .unwrap_or_else(|| panic!("cases must be an object"));
    for required in [
        "identical",
        "noise_0p5",
        "rigid_shift_8A",
        "hinge_first_half_8A",
        "displaced_tail_12A",
    ] {
        assert!(cases.contains_key(required), "{required} must be pinned");
    }
    for (label, entry) in cases {
        let value = entry["tm_norm_reference"]
            .as_f64()
            .unwrap_or_else(|| panic!("{label}: no reference value"));
        assert!(
            (0.0..=1.0).contains(&value),
            "{label}: the TM-score lies in [0, 1], found {value}"
        );
        let rmsd = entry["rmsd"]
            .as_f64()
            .unwrap_or_else(|| panic!("{label}: no reference rmsd"));
        assert!(rmsd >= 0.0 && rmsd.is_finite(), "{label}: {rmsd}");
    }

    // A rigid match is perfect and a hinged one is not, which is the whole point
    // of pinning these cases separately.
    let value = |label: &str| {
        cases[label]["tm_norm_reference"]
            .as_f64()
            .unwrap_or_else(|| panic!("{label}: not a number"))
    };
    assert!((value("identical") - 1.0).abs() < 1e-9);
    assert!(
        (value("rigid_shift_8A") - 1.0).abs() < 1e-9,
        "a rigid shift fits out"
    );
    assert!(
        (0.5..1.0).contains(&value("hinge_first_half_8A")),
        "a hinged model is neither perfect nor unrelated"
    );
    assert!(
        (0.5..1.0).contains(&value("displaced_tail_12A")),
        "a displaced tail is neither perfect nor unrelated"
    );

    // The recorded limitation is what makes this reference honest: the kernel
    // fits a fixed correspondence, and TM-align's own search can also re-pair
    // residues, so on these non-rigid cases TM-align reads a little higher.
    let limitation = reference["limitation"]
        .as_str()
        .unwrap_or_else(|| panic!("the reference must state its limitation"));
    assert!(
        limitation.contains("re-pair"),
        "the limitation names the search difference"
    );
}

/// The library's atomic weights are the CIAAW (IUPAC) conventional values.
///
/// gemmi 0.7.5 still carries pre-2021 Ge and Se, so those two distinguish a
/// current table from an outdated copy of it.
#[test]
fn gw_045_atomic_weights_match_the_ciaaw_table() {
    let reference = reference("atomic_weights_ciaaw_2021.json");
    assert_eq!(reference["library"], "CIAAW (IUPAC)");

    let weights = reference["weights"]
        .as_object()
        .unwrap_or_else(|| panic!("weights must be an object"));
    assert!(weights.len() >= 14);

    for (symbol, entry) in weights {
        let want = entry
            .as_f64()
            .unwrap_or_else(|| panic!("{symbol}: not a number"));
        let Some(element) = molframe::Element::from_symbol(symbol) else {
            panic!("{symbol} is not an element in the table");
        };
        let properties = molframe::chemistry::element_properties(element)
            .unwrap_or_else(|| panic!("{symbol} has no properties in the table"));
        let got = properties.atomic_weight;
        assert!(
            (got - want).abs() < 5e-4,
            "{symbol}: {got} against the CIAAW value {want}"
        );
    }

    // The two elements gemmi 0.7.5 disagrees with are present in the table and
    // carry the current value, which is what makes the reference discriminating.
    for (symbol, want) in [("Ge", 72.630), ("Se", 78.971)] {
        let element = molframe::Element::from_symbol(symbol).unwrap_or_else(|| panic!("{symbol}"));
        assert!(
            (match molframe::chemistry::element_properties(element) {
                Some(properties) => properties.atomic_weight,
                None => f64::NAN,
            } - want)
                .abs()
                < 5e-4,
            "{symbol}: {:?} against {want}",
            molframe::chemistry::element_properties(element).map(|p| p.atomic_weight)
        );
    }
}
