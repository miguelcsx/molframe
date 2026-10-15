use super::parse_lammps_dump;
use crate::FrameValue;

const ORTHOGONAL: &str = "ITEM: TIMESTEP
20
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0 10
0 20
0 30
ITEM: ATOMS id type x y z vx vy vz fx fy fz
2 1 2 3 4 0.2 0.3 0.4 2 3 4
1 1 1 2 3 0.1 0.2 0.3 1 2 3
";

#[test]
fn custom_columns_are_discovered_and_rows_are_sorted_by_id() {
    let frames =
        parse_lammps_dump(ORTHOGONAL).unwrap_or_else(|error| panic!("dump failed: {error}"));
    assert_eq!(frames.len(), 1);
    assert_eq!(frames[0].positions, [[1.0, 2.0, 3.0], [2.0, 3.0, 4.0]]);
    assert_eq!(
        frames[0].velocities,
        Some(vec![[0.1, 0.2, 0.3], [0.2, 0.3, 0.4]])
    );
    assert_eq!(
        frames[0].forces,
        Some(vec![[1.0, 2.0, 3.0], [2.0, 3.0, 4.0]])
    );
    assert_eq!(
        frames[0].data.get("lammps.step"),
        Some(&FrameValue::Integer(20))
    );
    let Some(cell) = frames[0].cell else {
        panic!("box should become a cell");
    };
    assert!(
        cell.lengths
            .iter()
            .zip([10.0, 20.0, 30.0])
            .all(|(observed, expected)| (observed - expected).abs() < f64::EPSILON)
    );
}

#[test]
fn scaled_restricted_triclinic_coordinates_are_cartesianised() {
    let source = "ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
1
ITEM: BOX BOUNDS xy xz yz pp pp pp
0 13 2
0 10 1
0 10 3
ITEM: ATOMS id xs ys zs
1 0.5 0.5 0.5
";
    let frames =
        parse_lammps_dump(source).unwrap_or_else(|error| panic!("triclinic dump failed: {error}"));
    let position = frames[0].positions[0];
    assert!((position[0] - 6.5).abs() < 1.0e-6);
    assert!((position[1] - 5.0).abs() < 1.0e-6);
    assert!((position[2] - 5.0).abs() < 1.0e-6);
    let Some(cell) = frames[0].cell else {
        panic!("box should become a cell");
    };
    assert!(cell.angles.iter().any(|angle| (*angle - 90.0).abs() > 1.0));
}

#[test]
fn changing_atom_count_is_refused_before_topology_can_drift() {
    let second = "ITEM: TIMESTEP
21
ITEM: NUMBER OF ATOMS
1
ITEM: BOX BOUNDS pp pp pp
0 10
0 20
0 30
ITEM: ATOMS id x y z
1 0 0 0
";
    assert!(parse_lammps_dump(&format!("{ORTHOGONAL}{second}")).is_err());
}
