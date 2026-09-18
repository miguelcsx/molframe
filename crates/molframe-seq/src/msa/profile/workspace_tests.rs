use super::*;

fn scoring() -> Scoring {
    Scoring {
        match_score: 3,
        mismatch_score: -2,
        gap_open: -5,
        gap_extend: -1,
    }
}

#[test]
fn packed_trace_round_trips_nibbles_at_exact_capacity() {
    let Ok(mut trace) = PackedTrace::new(31) else {
        panic!("small trace fits");
    };
    for (index, value) in (0_u8..16).cycle().take(31).enumerate() {
        trace.set(index, value);
    }
    assert_eq!(trace.values.len(), 16);
    for (index, value) in (0_u8..16).cycle().take(31).enumerate() {
        assert_eq!(trace.get(index), value);
    }
    assert_eq!(PackedTrace::required_bytes(100_000_000), Some(50_000_000));
}

#[test]
fn sparse_column_statistics_equal_pairwise_scoring() {
    let left = Profile {
        rows: vec![
            (0, b"AAC-".to_vec()),
            (1, b"ABC-".to_vec()),
            (2, b"-BCA".to_vec()),
        ],
    };
    let right = Profile {
        rows: vec![(3, b"AB-A".to_vec()), (4, b"BBCA".to_vec())],
    };
    let Ok(left_stats) = ColumnStats::new(&left) else {
        panic!("small profile fits");
    };
    let Ok(right_stats) = ColumnStats::new(&right) else {
        panic!("small profile fits");
    };
    for left_column in 0..left.width() {
        for right_column in 0..right.width() {
            let expected = left
                .rows
                .iter()
                .flat_map(|(_, left_row)| {
                    right.rows.iter().map(move |(_, right_row)| {
                        match (left_row[left_column], right_row[right_column]) {
                            (GAP, _) | (_, GAP) => 0_i64,
                            (first, second) => i64::from(scoring().substitution(first, second)),
                        }
                    })
                })
                .sum::<i64>();
            assert_eq!(
                left_stats.score(left_column, &right_stats, right_column, scoring()),
                Ok(expected)
            );
        }
    }
}

#[test]
fn rolling_workspace_preserves_affine_traceback_edges() {
    let left = Profile::singleton(0, b"ACGT");
    let right = Profile::singleton(1, b"ACXGT");
    let Ok(aligned) = align(&left, &right, scoring(), usize::MAX) else {
        panic!("small profiles align");
    };
    assert_eq!(aligned.rows[0].1, b"AC-GT");
    assert_eq!(aligned.rows[1].1, b"ACXGT");

    let empty = Profile::singleton(2, b"");
    let Ok(left_empty) = align(&empty, &left, scoring(), usize::MAX) else {
        panic!("empty left profile aligns");
    };
    assert_eq!(left_empty.rows[0].1, b"----");
    assert_eq!(left_empty.rows[1].1, b"ACGT");
    let Ok(right_empty) = align(&left, &empty, scoring(), usize::MAX) else {
        panic!("empty right profile aligns");
    };
    assert_eq!(right_empty.rows[0].1, b"ACGT");
    assert_eq!(right_empty.rows[1].1, b"----");
}
