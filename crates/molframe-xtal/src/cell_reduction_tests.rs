use super::*;

fn cell(lengths: [f64; 3], angles: [f64; 3]) -> UnitCell {
    UnitCell { lengths, angles }
}

/// The Gram matrix of a cell's edges.
fn metric(cell: &UnitCell) -> [[f64; 3]; 3] {
    let cosine = |degrees: f64| degrees.to_radians().cos();
    let [a, b, c] = cell.lengths;
    let [alpha, beta, gamma] = cell.angles;
    [
        [a * a, a * b * cosine(gamma), a * c * cosine(beta)],
        [a * b * cosine(gamma), b * b, b * c * cosine(alpha)],
        [a * c * cosine(beta), b * c * cosine(alpha), c * c],
    ]
}

fn determinant(matrix: [[f64; 3]; 3]) -> f64 {
    matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
}

fn integer_determinant(matrix: [[i32; 3]; 3]) -> i32 {
    matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
}

/// `Mᵀ G M`: the metric of the lattice whose edges are the columns of `basis`.
fn transformed(metric: [[f64; 3]; 3], basis: [[i32; 3]; 3]) -> [[f64; 3]; 3] {
    let mut result = [[0.0; 3]; 3];
    for (i, row) in result.iter_mut().enumerate() {
        for (j, entry) in row.iter_mut().enumerate() {
            for k in 0..3 {
                for l in 0..3 {
                    *entry += f64::from(basis[k][i]) * metric[k][l] * f64::from(basis[l][j]);
                }
            }
        }
    }
    result
}

fn assert_close(left: f64, right: f64, tolerance: f64) {
    assert!((left - right).abs() <= tolerance, "{left} != {right}");
}

/// A cell whose edges are `basis` columns of an orthorhombic cell.
fn skewed(basis: [[i32; 3]; 3]) -> UnitCell {
    let base = cell([10.0, 12.0, 15.0], [90.0; 3]);
    let gram = transformed(metric(&base), basis);
    let lengths = [gram[0][0].sqrt(), gram[1][1].sqrt(), gram[2][2].sqrt()];
    let angle = |first: usize, second: usize| {
        (gram[first][second] / (lengths[first] * lengths[second]))
            .acos()
            .to_degrees()
    };
    cell(lengths, [angle(1, 2), angle(0, 2), angle(0, 1)])
}

#[test]
fn a_reduced_cell_comes_back_unchanged() {
    let original = cell([10.0, 12.0, 15.0], [90.0; 3]);
    let reduced = niggli_reduce(&original, 1e-9, 100).expect("a valid cell");
    assert_eq!(reduced.change_of_basis, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
    assert!(reduced.converged);
    for (actual, expected) in reduced.lengths.into_iter().zip([10.0, 12.0, 15.0]) {
        assert_close(actual, expected, 1e-12);
    }
}

#[test]
fn a_skewed_setting_of_an_orthorhombic_lattice_is_recovered() {
    for basis in [
        [[1, 1, 0], [0, 1, 0], [0, 0, 1]],
        [[1, 0, 0], [3, 1, 0], [-2, 4, 1]],
        [[2, 3, 1], [1, 2, 1], [1, 1, 1]],
    ] {
        assert_eq!(integer_determinant(basis), 1, "{basis:?}");
        let reduced = niggli_reduce(&skewed(basis), 1e-9, 100).expect("a valid cell");
        assert!(reduced.converged, "{basis:?}");
        for (actual, expected) in reduced.lengths.into_iter().zip([10.0, 12.0, 15.0]) {
            assert_close(actual, expected, 1e-9);
        }
        for angle in reduced.angles {
            assert_close(angle, 90.0, 1e-7);
        }
    }
}

#[test]
fn the_change_of_basis_reproduces_the_reduced_metric() {
    let mut seed = 0x2545_f491_4f6c_dd1d_u64;
    let mut next = move || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        f64::from(u32::try_from(seed >> 40).expect("24 bits fit")) / f64::from(1_u32 << 24)
    };
    for _ in 0..200 {
        let original = cell(
            [
                5.0 + 40.0 * next(),
                5.0 + 40.0 * next(),
                5.0 + 40.0 * next(),
            ],
            [
                60.0 + 60.0 * next(),
                60.0 + 60.0 * next(),
                60.0 + 60.0 * next(),
            ],
        );
        let Ok(reduced) = niggli_reduce(&original, 1e-9, 100) else {
            continue;
        };
        assert!(reduced.converged, "{original:?}");
        assert_eq!(integer_determinant(reduced.change_of_basis), 1);
        let expected = transformed(metric(&original), reduced.change_of_basis);
        let actual = metric(&cell(reduced.lengths, reduced.angles));
        for (expected_row, actual_row) in expected.iter().zip(&actual) {
            for (left, right) in expected_row.iter().zip(actual_row) {
                assert_close(*left, *right, 1e-7 * left.abs().max(1.0));
            }
        }
        assert_close(
            determinant(actual),
            determinant(metric(&original)),
            1e-6 * determinant(actual),
        );
        assert!(reduced.lengths[0] <= reduced.lengths[1] + 1e-9);
        assert!(reduced.lengths[1] <= reduced.lengths[2] + 1e-9);
        let again = niggli_reduce(&cell(reduced.lengths, reduced.angles), 1e-9, 100)
            .expect("the reduced cell is valid");
        assert_eq!(again.change_of_basis, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
    }
}

#[test]
fn an_invalid_cell_or_tolerance_is_rejected() {
    let good = cell([10.0; 3], [90.0; 3]);
    assert_eq!(
        niggli_reduce(&good, -1.0, 10),
        Err(CellReductionError::InvalidTolerance)
    );
    assert_eq!(
        niggli_reduce(&good, f64::NAN, 10),
        Err(CellReductionError::InvalidTolerance)
    );
    assert_eq!(
        niggli_reduce(&cell([0.0, 1.0, 1.0], [90.0; 3]), 1e-9, 10),
        Err(CellReductionError::InvalidCell)
    );
    assert_eq!(
        niggli_reduce(&cell([1.0; 3], [10.0, 10.0, 170.0]), 1e-9, 10),
        Err(CellReductionError::InvalidCell)
    );
    assert_eq!(
        niggli_reduce(&cell([1.0, f64::NAN, 1.0], [90.0; 3]), 1e-9, 10),
        Err(CellReductionError::InvalidCell)
    );
}

#[test]
fn the_iteration_limit_stops_an_unfinished_reduction() {
    let reduced =
        niggli_reduce(&skewed([[1, 0, 0], [3, 1, 0], [-2, 4, 1]]), 1e-9, 1).expect("a valid cell");
    assert_eq!(reduced.iterations, 1);
    assert!(!reduced.converged);
}

/// A cell, its reduced parameters and its change of basis.
type Reference = ([f64; 3], [f64; 3], [f64; 6], [[i32; 3]; 3]);

/// Cells and results produced by Gemmi 0.7.5 (`GruberVector.niggli_reduce` on a
/// primitive lattice, default tolerance).
#[test]
fn reduced_cells_and_changes_of_basis_agree_with_the_gemmi_reference() {
    let cases: [Reference; 3] = [
        (
            [
                8.873_131_688_480_825,
                12.137_567_216_961_39,
                57.158_064_931_047_626,
            ],
            [
                99.750_687_423_710_63,
                79.519_449_898_383_27,
                90.911_201_744_261,
            ],
            [
                8.873_131_688_480_825,
                12.137_567_216_961_39,
                55.470_893_570_193_645,
                92.684_770_811_684_49,
                91.374_917_401_103_76,
                90.911_201_744_261,
            ],
            [[-1, 0, -1], [0, -1, 1], [0, 0, 1]],
        ),
        (
            [
                42.034_929_650_175_094,
                16.121_881_528_342_097,
                54.578_709_327_849_73,
            ],
            [
                67.371_860_593_482_25,
                52.645_974_990_194_3,
                66.061_516_332_897_54,
            ],
            [
                16.121_881_528_342_097,
                38.430_505_549_774_5,
                44.292_194_345_564_575,
                104.929_139_580_866_42,
                95.107_606_089_078_49,
                91.392_669_074_470_73,
            ],
            [[0, 1, -1], [-1, -1, 0], [0, 0, 1]],
        ),
        (
            [
                59.934_136_819_782_1,
                12.920_844_669_652_055,
                18.437_294_356_556_563,
            ],
            [
                78.577_577_022_362_16,
                54.870_936_232_735_03,
                119.630_793_363_466_7,
            ],
            [
                12.920_844_669_652_055,
                18.437_294_356_556_563,
                32.366_603_293_251_51,
                80.594_641_172_985_46,
                86.762_554_281_742_95,
                78.577_577_022_362_16,
            ],
            [[0, 0, 1], [1, 0, 3], [0, 1, -2]],
        ),
    ];
    for (lengths, angles, expected, basis) in cases {
        let reduced = niggli_reduce(&cell(lengths, angles), 1e-9, 100).expect("a valid cell");
        assert!(reduced.converged);
        assert_eq!(reduced.change_of_basis, basis, "{lengths:?} {angles:?}");
        let actual = [
            reduced.lengths[0],
            reduced.lengths[1],
            reduced.lengths[2],
            reduced.angles[0],
            reduced.angles[1],
            reduced.angles[2],
        ];
        for (left, right) in actual.into_iter().zip(expected) {
            assert_close(left, right, 1e-10);
        }
    }
}
