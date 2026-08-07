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
            if left_atom == right_atom
                || !finite(positions[left_atom as usize])
                || !finite(positions[right_atom as usize])
            {
                continue;
            }
            let squared = distance_squared(
                positions[left_atom as usize],
                positions[right_atom as usize],
                None,
            );
            if squared <= 2.25 {
                expected.push(NeighborPair::new(left_atom, right_atom, squared));
            }
        }
    }
    canonicalise(&mut expected);
    assert_eq!(actual, expected);
}
