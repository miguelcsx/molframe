use super::*;

#[test]
fn an_absent_band_visits_every_column() {
    assert_eq!(banded_columns(1, 10, None), (1, 10));
    assert_eq!(banded_columns(9, 10, None), (1, 10));
}

#[test]
fn a_band_restricts_the_columns_to_the_diagonal_neighbourhood() {
    assert_eq!(banded_columns(5, 100, Some(2)), (3, 8));
    assert_eq!(banded_columns(50, 100, Some(1)), (49, 52));
}

#[test]
fn a_band_is_clamped_to_the_row_at_both_ends() {
    assert_eq!(banded_columns(1, 100, Some(10)), (1, 12));
    assert_eq!(banded_columns(98, 100, Some(10)), (88, 100));
}

#[test]
fn the_visited_columns_are_exactly_those_the_band_admits() {
    // The band used to be tested inside a full traversal; this holds the
    // restricted range to the same set of cells that test would have kept.
    for width in [4_usize, 17, 64] {
        for band in [0_usize, 1, 3, 64] {
            for row in 1..width {
                let (first, last) = banded_columns(row, width, Some(band));
                let restricted: Vec<usize> = (first..last).collect();
                let filtered: Vec<usize> = (1..width)
                    .filter(|column| row.abs_diff(*column) <= band)
                    .collect();
                assert_eq!(restricted, filtered, "row {row} width {width} band {band}");
            }
        }
    }
}

#[test]
fn a_band_never_yields_a_reversed_range() {
    for row in 1..40_usize {
        for band in [0_usize, 1, 5] {
            let (first, last) = banded_columns(row, 20, Some(band));
            assert!(last >= first, "reversed at row {row} band {band}");
        }
    }
}

#[test]
fn narrow_band_stores_linear_cells_for_long_sequences() {
    let trace = layout(100_001, 100_001, Some(2)).expect("band");
    assert!(trace.bits.len() < 800_010);
    assert_eq!(trace.bits[0], ALL_STOP);
    assert_eq!(trace.at(50_000, 60_000), 0);
    let sequence = vec![b'A'; 100_000];
    let scoring = crate::Scoring {
        match_score: 2,
        mismatch_score: -1,
        gap_open: -2,
        gap_extend: -1,
    };
    let alignment = run(
        &sequence,
        &sequence,
        &scoring,
        -2,
        -1,
        Mode::Global,
        Some(2),
    )
    .expect("long band");
    assert_eq!(alignment.score, 200_000);
    assert_eq!(alignment.columns.len(), 100_000);
}

#[test]
fn wide_band_matches_dense_scores_and_traceback_for_all_modes() {
    let scoring = crate::Scoring {
        match_score: 2,
        mismatch_score: -1,
        gap_open: -2,
        gap_extend: -1,
    };
    for mode in [Mode::Global, Mode::Local, Mode::SemiGlobal] {
        for (left, right) in [
            (b"ACGTA".as_slice(), b"ACTA".as_slice()),
            (b"", b"AC"),
            (b"AC", b""),
        ] {
            let dense = run(left, right, &scoring, -2, -1, mode, None).expect("dense");
            let banded = run(left, right, &scoring, -2, -1, mode, Some(10)).expect("band");
            assert_eq!(dense, banded);
        }
    }
}

#[test]
fn banded_global_covers_both_sequences_after_a_leading_gap() {
    let scoring = crate::Scoring {
        match_score: 2,
        mismatch_score: -1,
        gap_open: -3,
        gap_extend: -1,
    };
    let left = b"AAAAB";
    let right = b"B";
    let banded = run(left, right, &scoring, -3, -1, Mode::Global, Some(4)).expect("band");
    let dense = crate::align::global(left, right, scoring).expect("dense");
    assert_eq!(banded.score, dense.score);
    let covered = |side: fn(&Column) -> Option<usize>| {
        banded.columns.iter().filter_map(side).collect::<Vec<_>>()
    };
    assert_eq!(covered(|c| c.left), vec![0, 1, 2, 3, 4]);
    assert_eq!(covered(|c| c.right), vec![0]);
}

#[test]
fn traceback_costs_one_byte_per_cell() {
    let trace = layout(1001, 2001, None).expect("dense");
    assert_eq!(trace.bits.len(), 1 + 1001 * 2001);
}
