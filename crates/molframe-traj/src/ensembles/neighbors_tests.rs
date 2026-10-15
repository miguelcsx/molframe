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

#[test]
fn borrowed_positions_match_the_timestep_adapter() {
    let frame = timestep(1.8);
    let mut from_timestep = FrameNeighborList::new([0], [1], 2.0, 1.0);
    let mut from_positions = FrameNeighborList::new([0], [1], 2.0, 1.0);

    let expected = from_timestep
        .pairs(&frame)
        .unwrap_or_else(|error| panic!("timestep query failed: {error}"));
    let actual = from_positions
        .pairs_positions(&frame.positions, frame.cell)
        .unwrap_or_else(|error| panic!("borrowed-coordinate query failed: {error}"));

    assert_eq!(actual, expected);
    assert_eq!(from_positions.statistics(), from_timestep.statistics());
}
