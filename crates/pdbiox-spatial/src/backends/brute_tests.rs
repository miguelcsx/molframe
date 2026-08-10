use super::*;

#[test]
fn simd_blocks_and_scalar_tail_match_the_direct_definition() {
    let positions = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 3.0, 0.0],
        [f32::NAN, 0.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let left = [0, 1, 2];
    let right = [0, 1, 2, 3, 4, 5];
    let actual = pairs(&positions, &left, &right, 2.25, None);
    let mut expected = Vec::new();
    for left_atom in left {
        for right_atom in right {
            let Ok(left_index) = usize::try_from(left_atom) else {
                panic!("test atom index must fit usize");
            };
            let Ok(right_index) = usize::try_from(right_atom) else {
                panic!("test atom index must fit usize");
            };
            if left_atom == right_atom
                || !finite(positions[left_index])
                || !finite(positions[right_index])
            {
                continue;
            }
            let squared = distance_squared(positions[left_index], positions[right_index], None);
            if squared <= 2.25 {
                expected.push(NeighborPair::new(left_atom, right_atom, squared));
            }
        }
    }
    canonicalise(&mut expected);
    assert_eq!(actual, expected);
}

#[test]
fn masked_simd_lanes_never_emit_pairs() {
    let positions = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [f32::NAN, 0.0, 0.0],
        [2.0, 0.0, 0.0],
    ];
    let actual = pairs(&positions, &[0], &[1, u32::MAX, 2, 3], 4.0, None);

    assert_eq!(
        actual,
        vec![NeighborPair::new(0, 1, 1.0), NeighborPair::new(0, 3, 4.0)]
    );
}
