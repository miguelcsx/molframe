use super::*;
use crate::space_group_by_hall;
use molframe_core::structure::UnitCell;

fn cubic(edge: f64) -> CellTransform {
    let cell = UnitCell {
        lengths: [edge, edge, edge],
        angles: [90.0; 3],
    };
    CellTransform::new(&cell).expect("a cubic cell is valid")
}

fn assert_close(left: f64, right: f64) {
    assert!(
        (left - right).abs() <= 1e-12 * left.abs().max(right.abs()).max(1.0),
        "{left} != {right}"
    );
}

#[test]
fn equal_steps_in_inverse_d2_put_every_limit_on_the_grid() {
    let spacings = [0.0, 1.0, 2.0, 3.0, 4.0];
    let binner = ResolutionBinner::new(BinMethod::Dstar2, 2, &spacings).expect("valid");
    // Four sub-shells of width one: limits 1, 2, 3 and infinity; shells keep
    // the odd ones.
    assert_eq!(binner.len(), 2);
    assert_close(binner.limits()[0], 2.0);
    assert!(binner.limits()[1].is_infinite());
    assert_close(binner.midpoints()[0], 1.0);
    assert_close(binner.midpoints()[1], 3.0);
    assert_eq!(binner.bin_index(2.0), 0);
    assert_eq!(binner.bin_index(2.0001), 1);
    assert_eq!(binner.bin_index(1.0e9), 1);
}

#[test]
fn every_method_orders_its_limits_and_ends_open() {
    let spacings: Vec<f64> = (1..=200).map(|n| f64::from(n) * 0.01).collect();
    for method in [
        BinMethod::EqualCount,
        BinMethod::Dstar,
        BinMethod::Dstar2,
        BinMethod::Dstar3,
    ] {
        let binner = ResolutionBinner::new(method, 8, &spacings).expect("valid");
        assert_eq!(binner.len(), 8);
        assert!(binner.limits().windows(2).all(|pair| pair[0] <= pair[1]));
        assert!(binner.limits()[7].is_infinite());
        assert!(binner.midpoints().windows(2).all(|pair| pair[0] <= pair[1]));
    }
}

#[test]
fn equal_count_shells_hold_nearly_equal_numbers_of_reflections() {
    let spacings: Vec<f64> = (0..1000).map(|n| (f64::from(n) / 1000.0).powi(2)).collect();
    let binner = ResolutionBinner::new(BinMethod::EqualCount, 10, &spacings).expect("valid");
    let mut counts = [0_usize; 10];
    for bin in binner.bin_indices(&spacings) {
        counts[bin] += 1;
    }
    assert!(
        counts.iter().all(|count| (95..=105).contains(count)),
        "{counts:?}"
    );
}

#[test]
fn shell_edges_report_resolution_in_angstrom() {
    let spacings = [0.01, 0.04, 0.09, 0.16];
    let binner = ResolutionBinner::new(BinMethod::Dstar2, 1, &spacings).expect("valid");
    assert_close(binner.d_max_of_bin(0).expect("bin 0"), 10.0);
    assert_close(binner.d_min_of_bin(0).expect("bin 0"), 2.5);
    assert!(binner.d_min_of_bin(1).is_none());
}

#[test]
fn invalid_binning_input_is_rejected() {
    let good = [0.1, 0.2];
    assert_eq!(
        ResolutionBinner::new(BinMethod::Dstar, 0, &good),
        Err(ReflectionBinningError::NoBins)
    );
    assert_eq!(
        ResolutionBinner::new(BinMethod::Dstar, 2, &[]),
        Err(ReflectionBinningError::NoReflections)
    );
    assert_eq!(
        ResolutionBinner::new(BinMethod::Dstar, 2, &[0.1, f64::NAN]),
        Err(ReflectionBinningError::InvalidSpacing)
    );
    assert_eq!(
        ResolutionBinner::new(BinMethod::Dstar, 2, &[-0.1]),
        Err(ReflectionBinningError::InvalidSpacing)
    );
}

fn lattice(limit: i32) -> Vec<[i32; 3]> {
    let mut out = Vec::new();
    for h in -limit..=limit {
        for k in -limit..=limit {
            for l in -limit..=limit {
                if [h, k, l] != [0, 0, 0] {
                    out.push([h, k, l]);
                }
            }
        }
    }
    out
}

#[test]
fn constant_amplitudes_normalize_to_unity_in_p1() {
    let cell = cubic(20.0);
    let p1 = space_group_by_hall("P 1").expect("P1").symmetry_set();
    let reflections = lattice(6);
    let amplitudes = vec![5.0; reflections.len()];
    let binner = ResolutionBinner::for_reflections(BinMethod::EqualCount, 6, &cell, &reflections)
        .expect("valid");
    let multipliers =
        amplitude_normalizers(&cell, &p1, &reflections, &amplitudes, &binner).expect("valid");
    for multiplier in multipliers {
        assert_close(multiplier * 5.0, 1.0);
    }
}

#[test]
fn multiplicity_lowers_the_multiplier_of_special_reflections() {
    let cell = UnitCell {
        lengths: [20.0, 22.0, 24.0],
        angles: [90.0; 3],
    };
    let cell = CellTransform::new(&cell).expect("valid");
    let p222 = space_group_by_hall("P 2 2").expect("P222").symmetry_set();
    // Epsilon is 1 for a general reflection and 2 on a twofold axis.
    let reflections = vec![[1, 2, 3], [0, 0, 3]];
    let amplitudes = vec![4.0; 2];
    let binner = ResolutionBinner::new(BinMethod::Dstar2, 1, &[0.01, 0.2]).expect("valid");
    let multipliers =
        amplitude_normalizers(&cell, &p222, &reflections, &amplitudes, &binner).expect("valid");
    assert_close(multipliers[0] / multipliers[1], 2.0_f64.sqrt());
}

#[test]
fn missing_amplitudes_stay_nan_and_do_not_shift_the_scale() {
    let cell = cubic(20.0);
    let p1 = space_group_by_hall("P 1").expect("P1").symmetry_set();
    let reflections = lattice(3);
    let mut amplitudes = vec![2.0; reflections.len()];
    amplitudes[0] = f64::NAN;
    let binner =
        ResolutionBinner::for_reflections(BinMethod::Dstar, 3, &cell, &reflections).expect("valid");
    let multipliers =
        amplitude_normalizers(&cell, &p1, &reflections, &amplitudes, &binner).expect("valid");
    assert!(multipliers[0].is_nan());
    assert_close(multipliers[1] * 2.0, 1.0);
}

#[test]
fn mismatched_columns_and_empty_symmetry_are_rejected() {
    let cell = cubic(20.0);
    let p1 = space_group_by_hall("P 1").expect("P1").symmetry_set();
    let binner = ResolutionBinner::new(BinMethod::Dstar, 1, &[0.1]).expect("valid");
    assert_eq!(
        amplitude_normalizers(&cell, &p1, &[[1, 0, 0]], &[], &binner),
        Err(ReflectionBinningError::LengthMismatch)
    );
    assert_eq!(
        amplitude_normalizers(
            &cell,
            &SymmetrySet::default(),
            &[[1, 0, 0]],
            &[1.0],
            &binner
        ),
        Err(ReflectionBinningError::NoOperations)
    );
}

/// Reference values produced by Gemmi 0.7.5 (`Binner` and
/// `calculate_amplitude_normalizers`) for the reflections below in `P 1 21 1`.
#[test]
fn shells_and_normalizers_agree_with_the_gemmi_reference() {
    let cell = CellTransform::new(&UnitCell {
        lengths: [31.5, 40.2, 52.7],
        angles: [90.0, 101.3, 90.0],
    })
    .expect("valid");
    let group = space_group_by_hall("P 2yb").expect("P21").symmetry_set();
    let mut reflections: Vec<[i32; 3]> = Vec::new();
    for h in -3..=3 {
        for k in 0..=3 {
            for l in -3..=3 {
                if [h, k, l] != [0, 0, 0] {
                    reflections.push([h, k, l]);
                }
            }
        }
    }
    assert_eq!(reflections.len(), 195);
    let amplitudes: Vec<f64> = reflections
        .iter()
        .map(|[h, k, l]| {
            let mix: i32 = h * 7 + k * 3 + l * 5;
            10.0 + f64::from(mix.rem_euclid(11))
        })
        .collect();
    let cases = [
        (
            BinMethod::Dstar2,
            [
                0.005_426_096_427_983_605,
                0.010_477_753_053_288_555,
                0.015_529_409_678_593_505,
            ],
            [
                0.065_032_433_060_513_12,
                0.064_964_494_804_662_17,
                0.065_246_770_452_518_14,
                0.065_503_342_893_295_87,
                0.065_743_777_998_462_22,
                0.046_487_871_243_535_59,
                0.064_699_663_922_063_06,
            ],
        ),
        (
            BinMethod::Dstar3,
            [
                0.008_207_638_829_023_61,
                0.012_986_461_667_431_31,
                0.016_998_559_405_218_83,
            ],
            [
                0.065_528_185_679_804_24,
                0.065_753_872_593_345_93,
                0.065_404_661_209_166_94,
                0.065_168_736_116_347_06,
                0.065_166_738_958_426_7,
                0.046_079_843_025_317_1,
                0.067_030_738_546_860_74,
            ],
        ),
    ];
    for (method, limits, expected) in cases {
        let binner =
            ResolutionBinner::for_reflections(method, 4, &cell, &reflections).expect("valid");
        for (actual, reference) in binner.limits().iter().zip(limits) {
            assert!(
                (actual - reference).abs() <= 1e-12 * reference,
                "{method:?}"
            );
        }
        let multipliers = amplitude_normalizers(&cell, &group, &reflections, &amplitudes, &binner)
            .expect("valid");
        for (row, reference) in [0, 7, 20, 40, 60, 100, 194].into_iter().zip(expected) {
            let actual = multipliers[row];
            assert!(
                (actual - reference).abs() <= 1e-9 * reference,
                "{method:?} row {row}: {actual} != {reference}"
            );
        }
    }
}

#[test]
fn every_binning_failure_has_a_registered_diagnostic_code() {
    use molframe_core::{Code, Diagnostic};
    for (error, code) in [
        (ReflectionBinningError::NoBins, Code::E5101),
        (ReflectionBinningError::NoReflections, Code::E5103),
        (ReflectionBinningError::LengthMismatch, Code::E5102),
        (ReflectionBinningError::TooLarge, Code::E1903),
    ] {
        let diagnostic = Diagnostic::from(error);
        assert_eq!(diagnostic.code(), code);
        assert!(diagnostic.code().is_registered());
    }
}
