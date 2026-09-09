use super::{
    PathFrameMetric, PathSimilarity, PathSimilarityError, frame_distance, path_similarity,
    path_similarity_view,
};
use std::mem::size_of;

#[test]
fn identical_ordered_paths_have_zero_distance() {
    let path = vec![vec![[0.0, 0.0, 0.0]], vec![[1.0, 0.0, 0.0]]];
    let Ok(similarity) = path_similarity(&path, &path, PathFrameMetric::CartesianRmsd, 1_024)
    else {
        panic!("valid paths");
    };
    assert!(similarity.hausdorff_distance.abs() < 1.0e-12);
    assert!(similarity.discrete_frechet_distance.abs() < 1.0e-12);
}

#[test]
fn frechet_retains_path_order_while_hausdorff_does_not() {
    let first = vec![vec![[0.0, 0.0, 0.0]], vec![[2.0, 0.0, 0.0]]];
    let second = vec![vec![[2.0, 0.0, 0.0]], vec![[0.0, 0.0, 0.0]]];
    let Ok(similarity) = path_similarity(&first, &second, PathFrameMetric::CartesianRmsd, 1_024)
    else {
        panic!("valid paths");
    };
    assert!(similarity.hausdorff_distance.abs() < 1.0e-12);
    assert!((similarity.discrete_frechet_distance - 2.0).abs() < 1.0e-12);
}

#[test]
fn path_allocation_requires_an_explicit_sufficient_ceiling() {
    let path = vec![vec![[0.0, 0.0, 0.0]]];
    assert_eq!(
        path_similarity(&path, &path, PathFrameMetric::CartesianRmsd, 0),
        Err(PathSimilarityError::MemoryLimit)
    );
}

#[test]
fn borrowed_paths_match_owned_paths() {
    let first = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.1, 0.0, 0.0],
        [1.1, 0.0, 0.0],
        [0.1, 1.0, 0.0],
    ];
    let second = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.2, 0.0, 0.0],
        [1.2, 0.0, 0.0],
        [0.2, 1.0, 0.0],
    ];
    let first_view = crate::FrameView::new(&first, 2, 3).expect("valid first view");
    let second_view = crate::FrameView::new(&second, 2, 3).expect("valid second view");
    let owned_first = vec![first[..3].to_vec(), first[3..].to_vec()];
    let owned_second = vec![second[..3].to_vec(), second[3..].to_vec()];
    assert_eq!(
        path_similarity_view(
            first_view,
            second_view,
            PathFrameMetric::CartesianRmsd,
            1024,
        ),
        path_similarity(
            &owned_first,
            &owned_second,
            PathFrameMetric::CartesianRmsd,
            1024,
        )
    );
}

#[test]
fn streaming_matches_the_matrix_reference_for_all_small_scalar_paths() {
    let paths = scalar_paths(4);
    for first in &paths {
        for second in &paths {
            let Ok(expected) = matrix_reference(first, second, PathFrameMetric::CartesianRmsd)
            else {
                panic!("the scalar reference path must be valid");
            };
            let Ok(actual) =
                path_similarity(first, second, PathFrameMetric::CartesianRmsd, usize::MAX)
            else {
                panic!("the scalar streaming path must be valid");
            };
            assert_eq!(
                actual,
                expected,
                "matrix mismatch for path lengths {} and {}",
                first.len(),
                second.len()
            );
        }
    }
}

#[test]
fn streaming_matches_the_matrix_reference_for_fitted_frames() {
    let first = vec![
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        vec![
            [0.1, 0.0, 0.0],
            [1.2, 0.1, 0.0],
            [0.0, 0.9, 0.1],
            [0.1, 0.0, 1.1],
        ],
    ];
    let second = vec![
        vec![
            [3.0, -2.0, 1.0],
            [3.0, -1.0, 1.0],
            [2.0, -2.0, 1.0],
            [3.0, -2.0, 2.0],
        ],
        vec![
            [0.2, 0.0, 0.0],
            [1.1, 0.0, 0.1],
            [0.0, 1.2, 0.0],
            [0.0, 0.1, 0.9],
        ],
        vec![
            [-1.0, 2.0, 0.5],
            [0.0, 2.0, 0.5],
            [-1.0, 3.0, 0.5],
            [-1.0, 2.0, 1.5],
        ],
    ];
    let Ok(expected) = matrix_reference(&first, &second, PathFrameMetric::FittedRmsd) else {
        panic!("the fitted matrix reference must succeed");
    };
    let Ok(actual) = path_similarity(&first, &second, PathFrameMetric::FittedRmsd, usize::MAX)
    else {
        panic!("the fitted streaming comparison must succeed");
    };
    assert_eq!(actual, expected);
}

#[test]
fn workspace_ceiling_tracks_columns_instead_of_frame_pairs() {
    let first = linear_path(257, 0.0);
    let second = linear_path(37, 0.0625);
    let workspace_bytes = 2 * second.len() * size_of::<f64>();
    let matrix_bytes = first.len() * second.len() * size_of::<f64>();
    assert!(workspace_bytes < matrix_bytes);
    assert_eq!(
        path_similarity(
            &first,
            &second,
            PathFrameMetric::CartesianRmsd,
            workspace_bytes - 1,
        ),
        Err(PathSimilarityError::MemoryLimit)
    );

    let Ok(expected) = matrix_reference(&first, &second, PathFrameMetric::CartesianRmsd) else {
        panic!("the scale reference must succeed");
    };
    let Ok(actual) = path_similarity(
        &first,
        &second,
        PathFrameMetric::CartesianRmsd,
        workspace_bytes,
    ) else {
        panic!("the streaming workspace must fit its exact ceiling");
    };
    assert_eq!(actual, expected);
}

#[test]
fn invalid_paths_are_rejected_before_the_workspace_ceiling() {
    let invalid = vec![vec![[f32::NAN, 0.0, 0.0]]];
    let valid = vec![vec![[0.0, 0.0, 0.0]]];
    assert_eq!(
        path_similarity(&invalid, &valid, PathFrameMetric::CartesianRmsd, 0),
        Err(PathSimilarityError::InvalidPaths)
    );
}

#[test]
fn fitted_degeneracy_is_still_reported_during_streaming() {
    let path = vec![vec![[0.0, 0.0, 0.0]; 4]];
    assert_eq!(
        path_similarity(&path, &path, PathFrameMetric::FittedRmsd, 16),
        Err(PathSimilarityError::DegenerateFit)
    );
}

fn scalar_paths(maximum_length: usize) -> Vec<Vec<Vec<[f32; 3]>>> {
    const VALUES: [f32; 3] = [-1.0, 0.0, 2.0];
    let mut paths = Vec::new();
    for length in 1..=maximum_length {
        let path_count = (0..length).fold(1_usize, |count, _| count * VALUES.len());
        for encoding in 0..path_count {
            let mut remainder = encoding;
            let mut path = Vec::with_capacity(length);
            for _ in 0..length {
                let value = VALUES[remainder % VALUES.len()];
                remainder /= VALUES.len();
                path.push(vec![[value, value * 0.25, value * -0.5]]);
            }
            paths.push(path);
        }
    }
    paths
}

fn linear_path(frame_count: usize, offset: f32) -> Vec<Vec<[f32; 3]>> {
    let mut coordinate = offset;
    (0..frame_count)
        .map(|_| {
            let frame = vec![[coordinate, coordinate * -0.25, coordinate * 0.5]];
            coordinate += 0.125;
            frame
        })
        .collect()
}

fn matrix_reference(
    first: &[Vec<[f32; 3]>],
    second: &[Vec<[f32; 3]>],
    metric: PathFrameMetric,
) -> Result<PathSimilarity, PathSimilarityError> {
    let rows = first.len();
    let columns = second.len();
    let mut distances = Vec::with_capacity(rows * columns);
    for first_frame in first {
        for second_frame in second {
            distances.push(frame_distance(first_frame, second_frame, metric)?);
        }
    }
    let first_to_second = (0..rows)
        .map(|row| {
            (0..columns)
                .map(|column| distances[row * columns + column])
                .fold(f64::INFINITY, f64::min)
        })
        .fold(0.0, f64::max);
    let second_to_first = (0..columns)
        .map(|column| {
            (0..rows)
                .map(|row| distances[row * columns + column])
                .fold(f64::INFINITY, f64::min)
        })
        .fold(0.0, f64::max);

    let mut previous: Vec<f64> = vec![0.0; columns];
    let mut current: Vec<f64> = vec![0.0; columns];
    for row in 0..rows {
        for column in 0..columns {
            let distance = distances[row * columns + column];
            current[column] = match (row, column) {
                (0, 0) => distance,
                (0, _) => current[column - 1].max(distance),
                (_, 0) => previous[0].max(distance),
                _ => previous[column]
                    .min(previous[column - 1])
                    .min(current[column - 1])
                    .max(distance),
            };
        }
        std::mem::swap(&mut previous, &mut current);
    }

    Ok(PathSimilarity {
        hausdorff_distance: first_to_second.max(second_to_first),
        discrete_frechet_distance: previous[columns - 1],
    })
}
