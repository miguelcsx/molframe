use super::{PathFrameMetric, path_similarity, path_similarity_view};

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
    assert!(path_similarity(&path, &path, PathFrameMetric::CartesianRmsd, 0).is_err());
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
