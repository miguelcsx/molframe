use super::knn_between;
use crate::{KnnOptions, SpatialBackend};
use molframe_core::ExecutionContext;

#[test]
fn separate_queries_index_into_the_targets() {
    let targets = [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
    let queries = [[4.0, 0.0, 0.0]];
    let context = ExecutionContext::default();
    let result = knn_between(
        &targets,
        &queries,
        2,
        KnnOptions::for_queries(None, None, true),
        SpatialBackend::Auto,
        None,
        &context,
    )
    .expect("a valid query succeeds");
    let (indices, distances) = result.into_columns();
    assert_eq!(indices, [1, 0]);
    assert_eq!(distances, [1.0, 4.0]);
}
