use super::*;

#[test]
fn bulk_measurements_match_scalar_kernels_and_mark_degenerate_rows() {
    let a = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
    let b = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
    let c = [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0]];
    let d = [[1.0, 1.0, 0.0], [1.0, 0.0, 0.0]];
    let mut values = [0.0; 2];

    angles_into(&a, &b, &c, &mut values).expect("matching inputs");
    assert_eq!(
        values[0].to_bits(),
        crate::angle(a[0], b[0], c[0]).expect("defined").to_bits()
    );
    assert!(values[1].is_nan());

    torsions_into(&a, &b, &c, &d, &mut values).expect("matching inputs");
    assert_eq!(
        values[0].to_bits(),
        crate::dihedral(a[0], b[0], c[0], d[0])
            .expect("defined")
            .to_bits()
    );
    assert!(values[1].is_nan());
}

#[test]
fn bulk_distances_reject_mismatched_storage() {
    let error =
        distances_into(&[[0.0; 3]], &[], &mut [0.0]).expect_err("mismatched inputs are rejected");
    assert_eq!(error.expected(), 1);
    assert_eq!(error.actual(), 0);
}

#[test]
fn every_short_tail_is_bit_identical_to_the_scalar_measurements() {
    let points = (0..12_i16)
        .map(|index| {
            let value = f32::from(index) + 0.25;
            [value, value * value + 1.0, value.mul_add(-0.75, 2.0)]
        })
        .collect::<Vec<_>>();
    for rows in 0..=9 {
        let first = &points[..rows];
        let second = &points[1..=rows];
        let third = &points[2..=rows + 1];
        let fourth = &points[3..=rows + 2];
        let mut distances = vec![0.0; rows];
        let mut angles = vec![0.0; rows];
        let mut torsions = vec![0.0; rows];
        distances_into(first, second, &mut distances).expect("distance batch");
        angles_into(first, second, third, &mut angles).expect("angle batch");
        torsions_into(first, second, third, fourth, &mut torsions).expect("torsion batch");
        for row in 0..rows {
            assert_eq!(
                distances[row].to_bits(),
                crate::distance(first[row], second[row]).to_bits()
            );
            assert_eq!(
                angles[row].to_bits(),
                crate::angle(first[row], second[row], third[row])
                    .expect("defined angle")
                    .to_bits()
            );
            assert_eq!(
                torsions[row].to_bits(),
                crate::dihedral(first[row], second[row], third[row], fourth[row])
                    .expect("defined torsion")
                    .to_bits()
            );
        }
    }
}
