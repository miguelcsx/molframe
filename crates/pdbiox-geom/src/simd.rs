//! Architecture-width kernels whose lane order preserves scalar semantics.

use num_traits::AsPrimitive;

#[cfg(target_arch = "x86_64")]
use wide::f64x4 as F64xN;
#[cfg(target_arch = "x86_64")]
const LANES: usize = 4;
#[cfg(not(target_arch = "x86_64"))]
use wide::f64x2 as F64xN;
#[cfg(not(target_arch = "x86_64"))]
const LANES: usize = 2;

pub(crate) fn squared_deviation_sum(left: &[[f32; 3]], right: &[[f32; 3]]) -> f64 {
    let mut total = 0.0;
    let mut left_blocks = left.chunks_exact(LANES);
    let mut right_blocks = right.chunks_exact(LANES);
    for (left_block, right_block) in (&mut left_blocks).zip(&mut right_blocks) {
        let squared = squared_block(left_block, right_block);
        for value in squared {
            total += value;
        }
    }
    for (&left, &right) in left_blocks.remainder().iter().zip(right_blocks.remainder()) {
        let dx = f64::from(left[0]) - f64::from(right[0]);
        let dy = f64::from(left[1]) - f64::from(right[1]);
        let dz = f64::from(left[2]) - f64::from(right[2]);
        total += dx * dx + dy * dy + dz * dz;
    }
    total
}

pub(crate) fn distances_into(left: &[[f32; 3]], right: &[[f32; 3]], output: &mut [f64]) {
    let mut left_blocks = left.chunks_exact(LANES);
    let mut right_blocks = right.chunks_exact(LANES);
    let mut output_blocks = output.chunks_exact_mut(LANES);
    for ((left, right), output) in (&mut left_blocks)
        .zip(&mut right_blocks)
        .zip(&mut output_blocks)
    {
        for (slot, squared) in output.iter_mut().zip(squared_block(left, right)) {
            *slot = squared.sqrt();
        }
    }
    for ((left, right), output) in left_blocks
        .remainder()
        .iter()
        .zip(right_blocks.remainder())
        .zip(output_blocks.into_remainder())
    {
        *output = crate::measure::distance(*left, *right);
    }
}

pub(crate) fn transform_all(
    positions: &mut [[f32; 3]],
    rotation: [[f64; 3]; 3],
    translation: [f64; 3],
) {
    let mut blocks = positions.chunks_exact_mut(LANES);
    for block in &mut blocks {
        transform_block(block, rotation, translation);
    }
    for position in blocks.into_remainder() {
        *position = transform_one(*position, rotation, translation);
    }
}

pub(crate) fn sum_positions(positions: &[[f32; 3]]) -> [f64; 3] {
    let mut sums = [0.0; 3];
    let mut blocks = positions.chunks_exact(LANES);
    for block in &mut blocks {
        let coordinates = coordinate_block(block);
        for ((x, y), z) in coordinates[0]
            .iter()
            .zip(&coordinates[1])
            .zip(&coordinates[2])
        {
            sums[0] += x;
            sums[1] += y;
            sums[2] += z;
        }
    }
    for position in blocks.remainder() {
        sums[0] += f64::from(position[0]);
        sums[1] += f64::from(position[1]);
        sums[2] += f64::from(position[2]);
    }
    sums
}

pub(crate) fn weighted_position_sum(positions: &[[f32; 3]], weights: &[f64]) -> ([f64; 3], f64) {
    let mut sums = [0.0; 3];
    let mut total = 0.0;
    let paired = positions.len().min(weights.len());
    let mut position_blocks = positions[..paired].chunks_exact(LANES);
    let mut weight_blocks = weights[..paired].chunks_exact(LANES);
    for (positions, weights) in (&mut position_blocks).zip(&mut weight_blocks) {
        let coordinates = coordinate_vectors(positions);
        let weight = F64xN::from(core::array::from_fn(|lane| weights[lane]));
        let products = coordinates.map(|coordinate| (coordinate * weight).to_array());
        for lane in 0..LANES {
            sums[0] += products[0][lane];
            sums[1] += products[1][lane];
            sums[2] += products[2][lane];
            total += weights[lane];
        }
    }
    for (position, weight) in position_blocks
        .remainder()
        .iter()
        .zip(weight_blocks.remainder())
    {
        sums[0] += f64::from(position[0]) * weight;
        sums[1] += f64::from(position[1]) * weight;
        sums[2] += f64::from(position[2]) * weight;
        total += weight;
    }
    (sums, total)
}

pub(crate) fn rmsf_update(
    frame: &[[f32; 3]],
    mean: &mut [[f64; 3]],
    summed_square: &mut [f64],
    inverse_seen: f64,
) {
    let mut offset = 0;
    while offset + LANES <= frame.len() {
        let points = coordinate_vectors(&frame[offset..offset + LANES]);
        let old: [F64xN; 3] = core::array::from_fn(|axis| {
            F64xN::from(core::array::from_fn(|lane| mean[offset + lane][axis]))
        });
        let before: [F64xN; 3] = core::array::from_fn(|axis| points[axis] - old[axis]);
        let updated: [F64xN; 3] =
            core::array::from_fn(|axis| old[axis] + before[axis] * F64xN::splat(inverse_seen));
        let after: [F64xN; 3] = core::array::from_fn(|axis| points[axis] - updated[axis]);
        let addition =
            (before[0] * after[0] + before[1] * after[1] + before[2] * after[2]).to_array();
        let arrays = updated.map(F64xN::to_array);
        for lane in 0..LANES {
            mean[offset + lane] = [arrays[0][lane], arrays[1][lane], arrays[2][lane]];
            summed_square[offset + lane] += addition[lane];
        }
        offset += LANES;
    }
    for atom in offset..frame.len() {
        let point = frame[atom].map(f64::from);
        let before: [f64; 3] = core::array::from_fn(|axis| point[axis] - mean[atom][axis]);
        for axis in 0..3 {
            mean[atom][axis] += before[axis] * inverse_seen;
        }
        let after: [f64; 3] = core::array::from_fn(|axis| point[axis] - mean[atom][axis]);
        summed_square[atom] += before[0] * after[0] + before[1] * after[1] + before[2] * after[2];
    }
}

pub(crate) fn inertia_components(
    positions: &[[f32; 3]],
    masses: Option<&[f64]>,
    centre: [f64; 3],
) -> [f64; 6] {
    let mut result = [0.0; 6];
    let mut offset = 0;
    while offset + LANES <= positions.len() {
        let coordinates = coordinate_vectors(&positions[offset..offset + LANES]);
        let x = coordinates[0] - F64xN::splat(centre[0]);
        let y = coordinates[1] - F64xN::splat(centre[1]);
        let z = coordinates[2] - F64xN::splat(centre[2]);
        let weight = match masses {
            Some(values) => F64xN::from(core::array::from_fn(|lane| values[offset + lane])),
            None => F64xN::splat(1.0),
        };
        let values = [
            weight * (y * y + z * z),
            -(weight * x * y),
            -(weight * x * z),
            weight * (x * x + z * z),
            -(weight * y * z),
            weight * (x * x + y * y),
        ]
        .map(F64xN::to_array);
        for (lane, _) in values[0].iter().enumerate() {
            for component in 0..6 {
                result[component] += values[component][lane];
            }
        }
        offset += LANES;
    }
    for row in offset..positions.len() {
        let weight = match masses {
            Some(values) => values[row],
            None => 1.0,
        };
        let x = f64::from(positions[row][0]) - centre[0];
        let y = f64::from(positions[row][1]) - centre[1];
        let z = f64::from(positions[row][2]) - centre[2];
        let values = [
            weight * (y * y + z * z),
            -weight * x * y,
            -weight * x * z,
            weight * (x * x + z * z),
            -weight * y * z,
            weight * (x * x + y * y),
        ];
        for component in 0..6 {
            result[component] += values[component];
        }
    }
    result
}

pub(crate) fn paired_moments(
    mobile: &[[f32; 3]],
    reference: &[[f32; 3]],
    mobile_centre: [f64; 3],
    reference_centre: [f64; 3],
    covariance: &mut [[f64; 3]; 3],
    mobile_scatter: &mut [[f64; 3]; 3],
    reference_scatter: &mut [[f64; 3]; 3],
) {
    let mut offset = 0;
    while offset + LANES <= mobile.len() {
        let mobile_points = coordinate_vectors(&mobile[offset..offset + LANES]);
        let reference_points = coordinate_vectors(&reference[offset..offset + LANES]);
        let mobile_delta: [F64xN; 3] =
            core::array::from_fn(|axis| mobile_points[axis] - F64xN::splat(mobile_centre[axis]));
        let reference_delta: [F64xN; 3] = core::array::from_fn(|axis| {
            reference_points[axis] - F64xN::splat(reference_centre[axis])
        });
        let covariance_values: [[[f64; LANES]; 3]; 3] = core::array::from_fn(|row| {
            core::array::from_fn(|column| (mobile_delta[row] * reference_delta[column]).to_array())
        });
        let mobile_values: [[[f64; LANES]; 3]; 3] = core::array::from_fn(|row| {
            core::array::from_fn(|column| (mobile_delta[row] * mobile_delta[column]).to_array())
        });
        let reference_values: [[[f64; LANES]; 3]; 3] = core::array::from_fn(|row| {
            core::array::from_fn(|column| {
                (reference_delta[row] * reference_delta[column]).to_array()
            })
        });
        for lane in 0..LANES {
            for row in 0..3 {
                for column in 0..3 {
                    covariance[row][column] += covariance_values[row][column][lane];
                    mobile_scatter[row][column] += mobile_values[row][column][lane];
                    reference_scatter[row][column] += reference_values[row][column][lane];
                }
            }
        }
        offset += LANES;
    }
    for (left, right) in mobile[offset..].iter().zip(&reference[offset..]) {
        let mobile_delta: [f64; 3] =
            core::array::from_fn(|axis| f64::from(left[axis]) - mobile_centre[axis]);
        let reference_delta: [f64; 3] =
            core::array::from_fn(|axis| f64::from(right[axis]) - reference_centre[axis]);
        for row in 0..3 {
            for column in 0..3 {
                covariance[row][column] += mobile_delta[row] * reference_delta[column];
                mobile_scatter[row][column] += mobile_delta[row] * mobile_delta[column];
                reference_scatter[row][column] += reference_delta[row] * reference_delta[column];
            }
        }
    }
}

pub(crate) fn squared_radius_sum(
    positions: &[[f32; 3]],
    centre: [f64; 3],
    weights: Option<&[f64]>,
) -> f64 {
    let paired = weights.map_or(positions.len(), |values| positions.len().min(values.len()));
    let mut sum = 0.0;
    let mut blocks = positions[..paired].chunks_exact(LANES);
    let mut offset = 0;
    for block in &mut blocks {
        let coordinates = coordinate_vectors(block);
        let dx = coordinates[0] - F64xN::splat(centre[0]);
        let dy = coordinates[1] - F64xN::splat(centre[1]);
        let dz = coordinates[2] - F64xN::splat(centre[2]);
        let mut squared = dx * dx + dy * dy + dz * dz;
        if let Some(values) = weights {
            squared *= F64xN::from(core::array::from_fn(|lane| values[offset + lane]));
        }
        for value in squared.to_array() {
            sum += value;
        }
        offset += LANES;
    }
    for (tail, position) in blocks.remainder().iter().enumerate() {
        let dx = f64::from(position[0]) - centre[0];
        let dy = f64::from(position[1]) - centre[1];
        let dz = f64::from(position[2]) - centre[2];
        let squared = dx * dx + dy * dy + dz * dz;
        sum += match weights {
            Some(values) => squared * values[offset + tail],
            None => squared,
        };
    }
    sum
}

fn squared_block(left: &[[f32; 3]], right: &[[f32; 3]]) -> [f64; LANES] {
    let left_x = F64xN::from(core::array::from_fn(|lane| f64::from(left[lane][0])));
    let left_y = F64xN::from(core::array::from_fn(|lane| f64::from(left[lane][1])));
    let left_z = F64xN::from(core::array::from_fn(|lane| f64::from(left[lane][2])));
    let right_x = F64xN::from(core::array::from_fn(|lane| f64::from(right[lane][0])));
    let right_y = F64xN::from(core::array::from_fn(|lane| f64::from(right[lane][1])));
    let right_z = F64xN::from(core::array::from_fn(|lane| f64::from(right[lane][2])));
    let dx = left_x - right_x;
    let dy = left_y - right_y;
    let dz = left_z - right_z;
    (dx * dx + dy * dy + dz * dz).to_array()
}

fn coordinate_vectors(block: &[[f32; 3]]) -> [F64xN; 3] {
    [
        F64xN::from(core::array::from_fn(|lane| f64::from(block[lane][0]))),
        F64xN::from(core::array::from_fn(|lane| f64::from(block[lane][1]))),
        F64xN::from(core::array::from_fn(|lane| f64::from(block[lane][2]))),
    ]
}

fn coordinate_block(block: &[[f32; 3]]) -> [[f64; LANES]; 3] {
    coordinate_vectors(block).map(F64xN::to_array)
}

fn transform_block(block: &mut [[f32; 3]], rotation: [[f64; 3]; 3], translation: [f64; 3]) {
    let x = F64xN::from(core::array::from_fn(|lane| f64::from(block[lane][0])));
    let y = F64xN::from(core::array::from_fn(|lane| f64::from(block[lane][1])));
    let z = F64xN::from(core::array::from_fn(|lane| f64::from(block[lane][2])));
    let arrays = [
        (x * F64xN::splat(rotation[0][0])
            + y * F64xN::splat(rotation[0][1])
            + z * F64xN::splat(rotation[0][2])
            + F64xN::splat(translation[0]))
        .to_array(),
        (x * F64xN::splat(rotation[1][0])
            + y * F64xN::splat(rotation[1][1])
            + z * F64xN::splat(rotation[1][2])
            + F64xN::splat(translation[1]))
        .to_array(),
        (x * F64xN::splat(rotation[2][0])
            + y * F64xN::splat(rotation[2][1])
            + z * F64xN::splat(rotation[2][2])
            + F64xN::splat(translation[2]))
        .to_array(),
    ];
    for lane in 0..LANES {
        block[lane] = [
            arrays[0][lane].as_(),
            arrays[1][lane].as_(),
            arrays[2][lane].as_(),
        ];
    }
}

fn transform_one(position: [f32; 3], rotation: [[f64; 3]; 3], translation: [f64; 3]) -> [f32; 3] {
    let point = position.map(f64::from);
    core::array::from_fn(|row| {
        let value = rotation[row][0] * point[0]
            + rotation[row][1] * point[1]
            + rotation[row][2] * point[2]
            + translation[row];
        value.as_()
    })
}
