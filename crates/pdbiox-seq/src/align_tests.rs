use super::{
    AlignError, Alignment, AlignmentMode, Column, RegionAlignError, RegionError, RegionOptions,
    align_region, global, global_banded, local, semi_global,
};
use crate::Scoring;

/// Rebuilds the two gapped rows of an alignment for inspection.
fn render(left: &[u8], right: &[u8], columns: &[Column]) -> (String, String) {
    let mut top = String::new();
    let mut bottom = String::new();
    for column in columns {
        match column.left {
            Some(index) => top.push(left[index] as char),
            None => top.push('-'),
        }
        match column.right {
            Some(index) => bottom.push(right[index] as char),
            None => bottom.push('-'),
        }
    }
    (top, bottom)
}

fn exact(result: Result<Alignment, AlignError>) -> Alignment {
    match result {
        Ok(alignment) => alignment,
        Err(error) => panic!("exact alignment failed: {error}"),
    }
}

#[test]
fn identical_sequences_align_with_no_gaps() {
    let alignment = exact(global(b"ACGT", b"ACGT", Scoring::simple()));
    assert_eq!(alignment.score, 4);
    assert_eq!(alignment.columns.len(), 4);
    assert!(
        alignment
            .columns
            .iter()
            .all(|column| column.left.is_some() && column.right.is_some())
    );
}

#[test]
fn a_global_alignment_pays_one_gap_to_match_the_rest() {
    let alignment = exact(global(b"ACGT", b"AGT", Scoring::simple()));
    assert_eq!(alignment.score, 1);
    let (top, bottom) = render(b"ACGT", b"AGT", &alignment.columns);
    assert_eq!(top, "ACGT");
    assert_eq!(bottom, "A-GT");
}

#[test]
fn a_local_alignment_finds_the_shared_core() {
    let alignment = exact(local(b"AACGTT", b"CGT", Scoring::simple()));
    assert_eq!(alignment.score, 3);
    let (top, bottom) = render(b"AACGTT", b"CGT", &alignment.columns);
    assert_eq!(top, "CGT");
    assert_eq!(bottom, "CGT");
}

#[test]
fn a_local_alignment_of_unrelated_sequences_is_empty() {
    let alignment = exact(local(b"AAAA", b"TTTT", Scoring::simple()));
    assert_eq!(alignment.score, 0);
    assert!(alignment.columns.is_empty());
}

#[test]
fn a_semi_global_alignment_ignores_the_overhang() {
    let alignment = exact(semi_global(b"TTACGTTT", b"ACGT", Scoring::simple()));
    assert_eq!(alignment.score, 4);
    let (top, bottom) = render(b"TTACGTTT", b"ACGT", &alignment.columns);
    assert_eq!(top, "ACGT");
    assert_eq!(bottom, "ACGT");
}

#[test]
fn a_wide_band_matches_the_unbanded_global_alignment() {
    let left = b"ACGTACGT";
    let right = b"AGTACT";
    let banded = exact(global_banded(left, right, Scoring::simple(), left.len()));
    let full = exact(global(left, right, Scoring::simple()));
    assert_eq!(banded, full);
}

#[test]
fn a_band_wide_enough_for_the_gap_still_scores_it() {
    // "ACGT" vs "AGT" needs one gap, a length difference of 1, so a band of 1
    // is enough to reach the corner and recover the score of 1.
    let banded = exact(global_banded(b"ACGT", b"AGT", Scoring::simple(), 1));
    assert_eq!(banded.score, 1);
}

#[test]
fn a_zero_band_keeps_equal_length_sequences_on_the_diagonal() {
    let banded = exact(global_banded(b"ACGT", b"ACGT", Scoring::simple(), 0));
    assert_eq!(banded.score, 4);
}

#[test]
fn a_band_without_a_complete_path_is_reported() {
    assert_eq!(
        global_banded(b"ACGT", b"AGT", Scoring::simple(), 0),
        Err(AlignError::NoAlignmentPath)
    );
}

#[test]
fn alignment_is_deterministic() {
    let first = exact(global(b"ACGTACGT", b"AGTACT", Scoring::simple()));
    let second = exact(global(b"ACGTACGT", b"AGTACT", Scoring::simple()));
    assert_eq!(first, second);
}

#[test]
fn a_region_alignment_reports_source_indices() {
    let options = RegionOptions {
        left: 2..6,
        right: 1..5,
        mode: AlignmentMode::Global,
        band: Some(0),
    };
    let result = align_region(b"XXACGTYY", b"ZACGTQ", Scoring::simple(), &options);
    let Ok(alignment) = result else {
        panic!("valid alignment regions");
    };
    assert_eq!(
        alignment.columns.first().and_then(|column| column.left),
        Some(2)
    );
    assert_eq!(
        alignment.columns.last().and_then(|column| column.right),
        Some(4)
    );
}

#[test]
fn an_out_of_bounds_region_is_rejected() {
    let start = 2;
    let end = 1;
    let options = RegionOptions {
        left: 0..4,
        right: start..end,
        mode: AlignmentMode::Local,
        band: None,
    };
    assert_eq!(
        align_region(b"AAAA", b"AAAA", Scoring::simple(), &options),
        Err(RegionAlignError::Region(RegionError {
            left: false,
            length: 4
        }))
    );
}

#[test]
fn exact_score_overflow_is_reported() {
    let scoring = Scoring {
        match_score: i32::MAX,
        mismatch_score: 0,
        gap_open: -1,
        gap_extend: -1,
    };
    assert_eq!(
        global(b"AA", b"AA", scoring),
        Err(AlignError::NumericOverflow)
    );
    assert_eq!(
        local(b"AA", b"AA", scoring),
        Err(AlignError::NumericOverflow)
    );
}

#[test]
fn score_only_matches_traceback_with_empty_and_unequal_inputs() {
    let sequences: &[&[u8]] = &[b"", b"A", b"AC", b"GATTACA", b"TAGACCA"];
    for left in sequences {
        for right in sequences {
            let scoring = Scoring::simple();
            assert_eq!(
                crate::global_score(left, right, scoring).expect("score"),
                global(left, right, scoring).expect("alignment").score,
            );
        }
    }
}
