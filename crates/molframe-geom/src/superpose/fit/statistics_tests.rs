use super::*;

/// The mean of a position set, computed the obvious way.
fn mean(set: &[[f32; 3]]) -> [f64; 3] {
    let count = exact_len(set.len());
    let mut centre = [0.0f64; 3];
    for point in set {
        for (axis, value) in centre.iter_mut().enumerate() {
            *value += f64::from(point[axis]);
        }
    }
    for value in &mut centre {
        *value /= count;
    }
    centre
}

/// The co-moment of two centred sets, summed directly with no blocking.
fn co_moment(left: &[[f32; 3]], right: &[[f32; 3]]) -> [[f64; 3]; 3] {
    let left_centre = mean(left);
    let right_centre = mean(right);
    let mut moment = [[0.0f64; 3]; 3];

    for (first, second) in left.iter().zip(right) {
        for (row, entries) in moment.iter_mut().enumerate() {
            let deviation = f64::from(first[row]) - left_centre[row];
            for (column, entry) in entries.iter_mut().enumerate() {
                *entry += deviation * (f64::from(second[column]) - right_centre[column]);
            }
        }
    }
    moment
}

fn assert_close(actual: &[[f64; 3]; 3], expected: &[[f64; 3]; 3], label: &str) {
    for (row, entries) in expected.iter().enumerate() {
        for (column, entry) in entries.iter().enumerate() {
            let difference = actual[row][column] - entry;
            assert!(
                difference.abs() < 1.0e-9,
                "{label}[{row}][{column}] was {} not {entry}",
                actual[row][column]
            );
        }
    }
}

#[test]
fn centres_are_the_means_of_each_set() {
    let mobile = [[0.0f32, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 3.0, 0.0]];
    let reference = [[1.0f32, 1.0, 1.0], [3.0, 1.0, 1.0], [2.0, 4.0, 1.0]];
    let statistics = fit_statistics(&mobile, &reference);

    for (axis, expected) in mean(&mobile).into_iter().enumerate() {
        assert!((statistics.mobile_centre[axis] - expected).abs() < 1.0e-12);
    }
    for (axis, expected) in mean(&reference).into_iter().enumerate() {
        assert!((statistics.reference_centre[axis] - expected).abs() < 1.0e-12);
    }
}

#[test]
fn scatter_and_covariance_match_the_directly_summed_moments() {
    let mobile = [
        [1.0f32, 2.0, 3.0],
        [4.0, 6.0, 1.0],
        [-2.0, 0.5, 7.0],
        [3.0, -1.0, 2.0],
        [0.0, 5.0, -4.0],
    ];
    let reference = [
        [2.0f32, 1.0, 0.0],
        [1.0, 3.0, 1.0],
        [0.0, 1.0, 5.0],
        [2.0, 2.0, 2.0],
        [-1.0, 4.0, 3.0],
    ];
    let statistics = fit_statistics(&mobile, &reference);

    assert_close(
        &statistics.mobile_scatter,
        &co_moment(&mobile, &mobile),
        "mobile_scatter",
    );
    assert_close(
        &statistics.reference_scatter,
        &co_moment(&reference, &reference),
        "reference_scatter",
    );
    assert_close(
        &statistics.covariance,
        &co_moment(&mobile, &reference),
        "covariance",
    );
}

#[test]
fn blocking_does_not_change_the_result_across_the_block_boundary() {
    // Spans several blocks, so the merge path is what produces the answer, and
    // it is compared against a direct sum that never blocks.
    let mobile: Vec<[f32; 3]> = (0..200i16)
        .map(|index| {
            let value = f32::from(index) * 0.37;
            [value, value.mul_add(-0.5, 3.0), (value * 0.25).sin() * 9.0]
        })
        .collect();
    let reference: Vec<[f32; 3]> = (0..200i16)
        .map(|index| {
            let value = f32::from(index) * 0.41;
            [value.mul_add(0.75, -1.0), value, (value * 0.5).cos() * 4.0]
        })
        .collect();

    let statistics = fit_statistics(&mobile, &reference);

    assert_close(
        &statistics.covariance,
        &co_moment(&mobile, &reference),
        "covariance",
    );
    assert_close(
        &statistics.mobile_scatter,
        &co_moment(&mobile, &mobile),
        "mobile_scatter",
    );
    for (axis, expected) in mean(&mobile).into_iter().enumerate() {
        assert!((statistics.mobile_centre[axis] - expected).abs() < 1.0e-9);
    }
}

#[test]
fn a_set_of_exactly_one_block_agrees_with_one_of_two_blocks() {
    let points: Vec<[f32; 3]> = (0..=i16::try_from(BLOCK_POINTS).expect("block fits"))
        .map(|index| {
            let value = f32::from(index);
            [value, value * 2.0, value * -0.5]
        })
        .collect();

    let Some(one_block) = points.get(..BLOCK_POINTS) else {
        panic!("test slice must exist");
    };
    let single = fit_statistics(one_block, one_block);
    assert_close(
        &single.mobile_scatter,
        &co_moment(one_block, one_block),
        "one block",
    );

    let spanning = fit_statistics(&points, &points);
    assert_close(
        &spanning.mobile_scatter,
        &co_moment(&points, &points),
        "two blocks",
    );
}

#[test]
fn unequal_lengths_pair_only_what_both_sets_have() {
    let mobile = [[1.0f32, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let reference = [[2.0f32, 0.0, 0.0], [0.0, 2.0, 0.0]];
    let statistics = fit_statistics(&mobile, &reference);

    let Some(paired) = mobile.get(..2) else {
        panic!("test slice must exist");
    };
    for (axis, expected) in mean(paired).into_iter().enumerate() {
        assert!((statistics.mobile_centre[axis] - expected).abs() < 1.0e-12);
    }
}

#[test]
fn an_empty_pair_produces_zero_statistics() {
    let statistics = fit_statistics(&[], &[]);

    assert!(statistics.mobile_centre.iter().all(|value| *value == 0.0));
    assert_close(&statistics.covariance, &[[0.0; 3]; 3], "covariance");
    assert_close(&statistics.mobile_scatter, &[[0.0; 3]; 3], "mobile_scatter");
}
