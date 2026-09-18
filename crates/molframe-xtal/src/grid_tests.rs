use super::{GridError, canonical_value_index, read_cube, read_dx};

/// Grid values are exact in the fixtures, but comparing floats by equality is
/// a habit worth not forming; a tolerance states the intent either way.
fn assert_values(actual: &[f32], expected: &[f32]) {
    assert_eq!(actual.len(), expected.len(), "value count");
    for (index, (actual, expected)) in actual.iter().zip(expected).enumerate() {
        assert!(
            (actual - expected).abs() < 1e-6,
            "value {index} was {actual}, expected {expected}"
        );
    }
}

/// A 2x2x2 cube in bohr whose values encode their own `x*100 + y*10 + z`
/// position, so a transposition error is visible in the value itself.
const CUBE: &str = "density
generated for a test
    2    0.000000    0.000000    0.000000
    2    1.000000    0.000000    0.000000
    2    0.000000    2.000000    0.000000
    2    0.000000    0.000000    4.000000
    8    8.000000    0.100000    0.200000    0.300000
    1    1.000000    1.100000    1.200000    1.300000
  0.0 1.0
  10.0 11.0
  100.0 101.0
  110.0 111.0
";

const DX: &str = "# a test potential
object 1 class gridpositions counts 2 2 2
origin 1.0 2.0 3.0
delta 0.5 0.0 0.0
delta 0.0 0.5 0.0
delta 0.0 0.0 0.5
object 2 class gridconnections counts 2 2 2
object 3 class array type double rank 0 items 8 data follows
0.0 1.0 10.0
11.0 100.0 101.0
110.0 111.0
attribute \"dep\" string \"positions\"
";

#[test]
fn a_cube_is_stored_first_axis_fastest_after_reading() {
    let Ok(grid) = read_cube(CUBE.as_bytes()) else {
        panic!("the fixture cube reads")
    };
    assert_eq!(grid.map.dimensions, [2, 2, 2]);
    // Written z-fastest as 0,1,10,11,100,101,110,111; stored x-fastest the
    // first pair must be the two x neighbours at y=0, z=0.
    assert_values(&grid.map.values[..4], &[0.0, 100.0, 10.0, 110.0]);
}

#[test]
fn cube_geometry_is_converted_from_bohr_to_angstroms() {
    let Ok(grid) = read_cube(CUBE.as_bytes()) else {
        panic!("the fixture cube reads")
    };
    // Two voxels of one bohr along x spans two bohr.
    let expected = 2.0 * 0.529_177_210_903;
    assert!(
        (grid.map.cell.lengths[0] - expected).abs() < 1e-6,
        "x span was {}",
        grid.map.cell.lengths[0]
    );
    for angle in grid.map.cell.angles {
        assert!((angle - 90.0).abs() < 1e-4, "axes are orthogonal: {angle}");
    }
}

#[test]
fn cube_nuclei_are_read_and_converted() {
    let Ok(grid) = read_cube(CUBE.as_bytes()) else {
        panic!("the fixture cube reads")
    };
    assert_eq!(grid.atoms.len(), 2);
    assert_eq!(grid.atoms[0].number, 8);
    let expected = 0.1 * 0.529_177_210_903;
    assert!((grid.atoms[0].position[0] - expected).abs() < 1e-9);
    assert_eq!(grid.fields, 1);
    assert!(grid.orbitals.is_empty());
}

#[test]
fn a_cube_written_in_angstroms_is_not_converted_again() {
    let angstrom = CUBE.replace(
        "    2    1.000000    0.000000    0.000000",
        "   -2    1.000000    0.000000    0.000000",
    );
    let Ok(grid) = read_cube(angstrom.as_bytes()) else {
        panic!("a negative voxel count marks angstroms")
    };
    assert!((grid.map.cell.lengths[0] - 2.0).abs() < 1e-6);
}

#[test]
fn a_truncated_cube_is_refused_rather_than_zero_filled() {
    let short = CUBE.replace("  110.0 111.0\n", "");
    assert!(matches!(
        read_cube(short.as_bytes()),
        Err(GridError::Truncated { .. })
    ));
}

#[test]
fn a_cube_value_that_is_not_a_number_reports_where_it_stopped() {
    let broken = CUBE.replace("  100.0 101.0", "  100.0 wrong");
    assert!(matches!(
        read_cube(broken.as_bytes()),
        Err(GridError::InvalidNumber { index: 5, .. })
    ));
}

#[test]
fn a_dx_grid_is_stored_first_axis_fastest_after_reading() {
    let Ok(map) = read_dx(DX.as_bytes()) else {
        panic!("the fixture grid reads")
    };
    assert_eq!(map.dimensions, [2, 2, 2]);
    assert_values(&map.values[..4], &[0.0, 100.0, 10.0, 110.0]);
    for (actual, expected) in map.origin.iter().zip([1.0, 2.0, 3.0]) {
        assert!((actual - expected).abs() < 1e-9, "origin was {actual}");
    }
}

#[test]
fn dx_deltas_become_the_span_of_the_whole_grid() {
    let Ok(map) = read_dx(DX.as_bytes()) else {
        panic!("the fixture grid reads")
    };
    for length in map.cell.lengths {
        assert!((length - 1.0).abs() < 1e-6, "span was {length}");
    }
}

#[test]
fn a_dx_grid_without_three_deltas_is_refused() {
    let irregular = DX.replace("delta 0.0 0.0 0.5\n", "");
    assert!(matches!(
        read_dx(irregular.as_bytes()),
        Err(GridError::InvalidHeader { .. })
    ));
}

#[test]
fn a_dx_item_count_that_disagrees_with_the_grid_is_refused() {
    let inconsistent = DX.replace("items 8", "items 9");
    assert!(matches!(
        read_dx(inconsistent.as_bytes()),
        Err(GridError::InvalidHeader { .. })
    ));
}

#[test]
fn an_orbital_cube_keeps_every_field_and_can_extract_one() {
    let orbital = "orbitals
two of them
   -1    0.000000    0.000000    0.000000
    2    1.000000    0.000000    0.000000
    1    0.000000    1.000000    0.000000
    1    0.000000    0.000000    1.000000
    8    8.000000    0.000000    0.000000    0.000000
    2    3    5
  0.0 1.0
  2.0 3.0
";
    let Ok(grid) = read_cube(orbital.as_bytes()) else {
        panic!("the fixture orbital cube reads")
    };
    assert_eq!(grid.fields, 2);
    assert_eq!(grid.orbitals, vec![3, 5]);
    let Some(second) = grid.field(1) else {
        panic!("the second field exists")
    };
    assert_values(&second.values, &[1.0, 3.0]);
    assert!(grid.field(2).is_none());
}

#[test]
fn rectangular_multifield_indices_preserve_every_scalar() {
    let counts = [2, 3, 4];
    let fields = 2;
    let mut destinations = (0..48)
        .map(|source| canonical_value_index(source, counts, fields))
        .collect::<Vec<_>>();
    destinations.sort_unstable();
    assert_eq!(destinations, (0..48).collect::<Vec<_>>());
    assert_eq!(canonical_value_index(0, counts, fields), 0);
    assert_eq!(canonical_value_index(2, counts, fields), 12);
    assert_eq!(canonical_value_index(8, counts, fields), 4);
    assert_eq!(canonical_value_index(24, counts, fields), 2);
}
