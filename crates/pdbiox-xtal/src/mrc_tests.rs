use super::*;

fn example() -> DensityMap {
    DensityMap {
        dimensions: [2, 2, 2],
        starts: [0; 3],
        sampling: [2, 2, 2],
        cell: UnitCell {
            lengths: [2.0, 2.0, 2.0],
            angles: [90.0; 3],
        },
        origin: [0.0; 3],
        space_group: 1,
        labels: vec!["pdbiox test".into()],
        extended_header: b"x,y,z".to_vec(),
        values: vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
    }
}

#[test]
fn mode_two_round_trip_preserves_map_and_metadata() {
    let expected = example();
    let bytes = expected.to_mrc_bytes().expect("map should encode");
    let actual = DensityMap::from_mrc_bytes(&bytes).expect("map should decode");
    assert_eq!(actual, expected);
}

#[test]
fn trilinear_sampling_uses_cartesian_cell_coordinates() {
    let map = example();
    let value = map
        .sample_cartesian([0.5, 0.5, 0.5], MapBoundary::Missing)
        .expect("coordinate is inside map");
    assert!((value - 3.5).abs() < f32::EPSILON);
    assert!(
        map.sample_cartesian([-1.0, 0.0, 0.0], MapBoundary::Missing)
            .is_none()
    );
}

#[test]
fn periodic_sampling_wraps_upper_neighbors() {
    let map = example();
    let value = map
        .sample_grid([1.5, 0.0, 0.0], MapBoundary::Periodic)
        .expect("periodic point should wrap");
    assert!((value - 0.5).abs() < f32::EPSILON);
}

#[test]
fn cubic_sampling_reconstructs_values_at_grid_points() {
    let map = example();
    let value = map
        .sample_grid_cubic([1.0, 1.0, 1.0], MapBoundary::Missing)
        .expect("exact grid point needs no outside samples");
    assert!((value - 7.0).abs() < f32::EPSILON);
}

#[test]
fn permuted_file_axes_are_canonicalised() {
    let map = example();
    let mut bytes = map.to_mrc_bytes().expect("map should encode");
    bytes[64..68].copy_from_slice(&2_i32.to_le_bytes());
    bytes[68..72].copy_from_slice(&1_i32.to_le_bytes());
    bytes[72..76].copy_from_slice(&3_i32.to_le_bytes());
    let actual = DensityMap::from_mrc_bytes(&bytes).expect("permuted map should decode");
    assert_eq!(actual.value([0, 1, 0]), Some(1.0));
    assert_eq!(actual.value([1, 0, 0]), Some(2.0));
}

#[test]
fn complex_transform_modes_are_refused() {
    let map = example();
    let mut bytes = map.to_mrc_bytes().expect("map should encode");
    bytes[12..16].copy_from_slice(&4_i32.to_le_bytes());
    assert_eq!(
        DensityMap::from_mrc_bytes(&bytes),
        Err(MrcError::UnsupportedMode(4))
    );
}

#[test]
fn external_ccpem_permuted_fixture_when_configured() {
    let Some(path) = std::env::var_os("PDBIOX_MRC_FIXTURE") else {
        return;
    };
    let bytes = std::fs::read(path).expect("configured CCP-EM fixture should be readable");
    let map = DensityMap::from_mrc_bytes(&bytes).expect("CCP-EM fixture should parse");
    assert_eq!(map.dimensions, [43, 25, 73]);
    assert_eq!(map.extended_header.len(), 160);
    let first = map.value([0, 0, 0]).expect("first density should exist");
    let second = map.value([0, 0, 1]).expect("second density should exist");
    assert!((first - 0.042_834_47).abs() < 1.0e-7);
    assert!((second - 0.026_947_163).abs() < 1.0e-7);
}
