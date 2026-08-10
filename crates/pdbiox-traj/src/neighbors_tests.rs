use super::FrameNeighborList;
use crate::Timestep;

fn timestep(distance: f32) -> Timestep {
    Timestep {
        positions: vec![[0.0, 0.0, 0.0], [distance, 0.0, 0.0]],
        ..Timestep::default()
    }
}

#[test]
fn candidates_are_reused_until_half_the_skin_is_exhausted() {
    let mut neighbors = FrameNeighborList::new([0], [1], 2.0, 1.0);
    assert_eq!(
        neighbors
            .pairs(&timestep(1.8))
            .unwrap_or_else(|error| panic!("first frame failed: {error}"))
            .len(),
        1
    );
    assert_eq!(
        neighbors
            .pairs(&timestep(2.1))
            .unwrap_or_else(|error| panic!("second frame failed: {error}"))
            .len(),
        0
    );
    assert_eq!(neighbors.statistics().rebuilds, 1);

    let _ = neighbors
        .pairs(&timestep(3.0))
        .unwrap_or_else(|error| panic!("third frame failed: {error}"));
    assert_eq!(neighbors.statistics().rebuilds, 2);
    assert_eq!(neighbors.statistics().frames, 3);
}
