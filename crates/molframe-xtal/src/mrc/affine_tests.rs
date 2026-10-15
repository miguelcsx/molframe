use super::*;
use crate::{MapBoundary, MrcBlockOptions, MrcBlockReader};
use std::io::Cursor;

fn map(angles: [f64; 3], origin: [f64; 3]) -> DensityMap {
    DensityMap {
        dimensions: [3, 3, 3],
        starts: [-2, 4, 1],
        sampling: [5, 6, 7],
        cell: UnitCell {
            lengths: [20.0, 30.0, 40.0],
            angles,
        },
        origin,
        space_group: 1,
        labels: Vec::new(),
        extended_header: Vec::new(),
        values: (0_u16..27).map(f32::from).collect(),
    }
}

fn apply(matrix: [[f64; 4]; 4], voxel: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|axis| {
        matrix[axis][3]
            + (0..3)
                .map(|column| matrix[axis][column] * voxel[column])
                .sum::<f64>()
    })
}

fn verify(bytes: Vec<u8>) {
    let map = DensityMap::from_mrc_bytes(&bytes).expect("map parses");
    let reader =
        MrcBlockReader::new(Cursor::new(bytes), MrcBlockOptions::default()).expect("header parses");
    let matrix = reader.descriptor().voxel_to_world().expect("valid affine");
    assert_eq!(matrix, map.voxel_to_world().expect("valid map affine"));
    let transform = CellTransform::new(&map.cell).expect("valid cell");
    let sampler = map.sampler().expect("valid sampler");
    for z in 0_u16..3 {
        for y in 0_u16..3 {
            for x in 0_u16..3 {
                let voxel = [x, y, z].map(f64::from);
                let world = apply(matrix, voxel);
                let has_origin = map.origin.iter().any(|value| value.abs() > f64::EPSILON);
                let expected = transform.to_cartesian(std::array::from_fn(|axis| {
                    (voxel[axis]
                        + if has_origin {
                            0.0
                        } else {
                            f64::from(map.starts[axis])
                        })
                        / usize_to_f64(map.sampling[axis])
                }));
                for axis in 0..3 {
                    let offset = if has_origin { map.origin[axis] } else { 0.0 };
                    assert!((world[axis] - expected[axis] - offset).abs() < 1e-9);
                }
                let stored = map.value([x, y, z].map(usize::from)).expect("stored voxel");
                // Interior samples avoid outside-domain rounding at exact edges.
                if [x, y, z] == [1, 1, 1] {
                    assert!(
                        (sampler
                            .sample_cartesian(world, MapBoundary::Missing)
                            .expect("inside")
                            - stored)
                            .abs()
                            < 1e-5
                    );
                }
            }
        }
    }
}

#[test]
fn orthogonal_and_skewed_affines_agree_with_sampling_and_start_indices() {
    for angles in [[90.0; 3], [75.0, 82.0, 68.0]] {
        verify(map(angles, [0.0; 3]).to_mrc_bytes().expect("encodes"));
    }
}

#[test]
fn explicit_origin_replaces_starts_even_with_permuted_file_axes() {
    let mut bytes = map([75.0, 82.0, 68.0], [1.25, -3.5, 7.0])
        .to_mrc_bytes()
        .expect("encodes");
    bytes[64..68].copy_from_slice(&2_i32.to_le_bytes());
    bytes[68..72].copy_from_slice(&3_i32.to_le_bytes());
    bytes[72..76].copy_from_slice(&1_i32.to_le_bytes());
    verify(bytes);
}

#[test]
fn invalid_caller_constructed_geometry_is_rejected() {
    let mut map = map([90.0; 3], [0.0; 3]);
    map.sampling[0] = 0;
    assert_eq!(map.voxel_to_world(), Err(MrcError::InvalidHeader));
    map.sampling[0] = 5;
    map.origin[0] = f64::NAN;
    assert_eq!(map.voxel_to_world(), Err(MrcError::InvalidHeader));
    map.origin[0] = 0.0;
    map.cell.angles[0] = 0.0;
    assert_eq!(map.voxel_to_world(), Err(MrcError::InvalidHeader));
}
