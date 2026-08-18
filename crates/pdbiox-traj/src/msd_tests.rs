use super::{mean_squared_displacement, mean_squared_displacement_view};

#[test]
fn windowed_msd_averages_atoms_and_time_origins() {
    let frames = vec![
        vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        vec![[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]],
        vec![[2.0, 0.0, 0.0], [0.0, 4.0, 0.0]],
    ];
    let Ok(msd) = mean_squared_displacement(&frames, &[], 2) else {
        panic!("valid trajectory");
    };
    assert!(msd[0].value.abs() < 1.0e-12);
    assert!((msd[1].value - 2.5).abs() < 1.0e-12);
    assert!((msd[2].value - 10.0).abs() < 1.0e-12);
    assert_eq!(msd[1].observations, 4);
}

#[test]
fn selection_and_unwrapping_are_explicit() {
    let frames = vec![
        vec![[0.0; 3], [0.0; 3]],
        vec![[3.0, 0.0, 0.0], [9.0, 0.0, 0.0]],
    ];
    let Ok(msd) = mean_squared_displacement(&frames, &[0], 1) else {
        panic!("valid selection");
    };
    assert!((msd[1].value - 9.0).abs() < 1.0e-12);
}

#[test]
fn borrowed_frame_views_match_owned_msd_without_an_atom_index_copy() {
    let frames = vec![
        vec![[0.0_f32, 0.0, 0.0], [0.0, 0.0, 0.0]],
        vec![[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]],
        vec![[2.0, 0.0, 0.0], [0.0, 4.0, 0.0]],
    ];
    let positions = frames.iter().flatten().copied().collect::<Vec<[f32; 3]>>();
    let view = crate::FrameView::new(&positions, 3, 2).expect("contiguous test coordinates");
    assert_eq!(
        mean_squared_displacement_view(view, &[], 2),
        mean_squared_displacement(&frames, &[], 2)
    );
}
