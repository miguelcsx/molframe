use super::{PathFrameMetric, path_similarity};

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
