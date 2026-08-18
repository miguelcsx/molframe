use super::{FrameView, FrameViewError};

#[test]
fn frame_view_borrows_contiguous_frames_without_repacking() {
    let positions = [
        [0.0_f32, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
    ];
    let view = FrameView::new(&positions, 2, 2).expect("matching dimensions");
    let Some(first) = view.frame(0) else {
        panic!("first frame is present");
    };
    let Some(second) = view.frame(1) else {
        panic!("second frame is present");
    };
    assert_eq!(first.as_ptr(), positions.as_ptr());
    assert_eq!(second.as_ptr(), positions[2..].as_ptr());
    assert_eq!(second, &positions[2..]);
}

#[test]
fn frame_view_rejects_inconsistent_or_overflowing_dimensions() {
    let positions = [[0.0_f32, 0.0, 0.0]];
    assert_eq!(
        FrameView::new(&positions, 1, 2).err(),
        Some(FrameViewError::DimensionMismatch)
    );
    assert_eq!(
        FrameView::new(&positions, usize::MAX, 2).err(),
        Some(FrameViewError::DimensionOverflow)
    );
}
