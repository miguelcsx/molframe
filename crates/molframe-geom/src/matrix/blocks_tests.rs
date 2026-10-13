use super::*;
use molframe_core::ExecutionContext;

#[test]
fn blocks_cover_every_row_once_with_a_short_tail() {
    let blocks: Vec<_> = RowBlocks::new(7, 3).expect("positive rows").collect();
    assert_eq!(blocks, [0..3, 3..6, 6..7]);
    assert!(RowBlocks::new(7, 0).is_none());
    assert_eq!(RowBlocks::new(0, 4).expect("positive rows").count(), 0);
}

#[test]
fn a_block_matches_the_matching_rows_of_the_full_matrix() {
    let points = [[0.0, 0.0, 0.0], [3.0, 4.0, 0.0], [0.0, 0.0, 12.0]];
    let context = ExecutionContext::default();
    let block = distance_rows_with_context(&points, 1..3, &context).expect("fits");
    assert_eq!(block.len(), 6);
    assert!((block[0] - 5.0).abs() < 1e-12);
    assert!((block[3] - 12.0).abs() < 1e-12);
    let narrow = distance_rows_f32_with_context(&points, 0..1, &context).expect("fits");
    assert!(
        narrow
            .iter()
            .zip([0.0, 5.0, 12.0])
            .all(|(a, b)| (a - b).abs() < 1e-6)
    );
}
