//! Vector arithmetic for caller-owned angle and torsion batches.

#[cfg(target_arch = "x86_64")]
use wide::f64x4 as F64xN;
#[cfg(target_arch = "x86_64")]
const LANES: usize = 4;
#[cfg(not(target_arch = "x86_64"))]
use wide::f64x2 as F64xN;
#[cfg(not(target_arch = "x86_64"))]
const LANES: usize = 2;

pub(crate) fn angles_into(
    first: &[[f32; 3]],
    vertices: &[[f32; 3]],
    third: &[[f32; 3]],
    output: &mut [f64],
) {
    let mut offset = 0;
    while offset + LANES <= first.len() {
        let a = displacement_vectors(&vertices[offset..], &first[offset..]);
        let b = displacement_vectors(&vertices[offset..], &third[offset..]);
        let aa = dot(a, a).to_array();
        let bb = dot(b, b).to_array();
        let ab = dot(a, b).to_array();
        for lane in 0..LANES {
            output[offset + lane] = if aa[lane] <= 0.0 || bb[lane] <= 0.0 {
                f64::NAN
            } else {
                (ab[lane] * (aa[lane] * bb[lane]).sqrt().recip())
                    .clamp(-1.0, 1.0)
                    .acos()
            };
        }
        offset += LANES;
    }
    for row in offset..first.len() {
        output[row] = match crate::measure::angle(first[row], vertices[row], third[row]) {
            Some(angle) => angle,
            None => f64::NAN,
        };
    }
}

pub(crate) fn torsions_into(
    first: &[[f32; 3]],
    second: &[[f32; 3]],
    third: &[[f32; 3]],
    fourth: &[[f32; 3]],
    output: &mut [f64],
) {
    let mut offset = 0;
    while offset + LANES <= first.len() {
        let a = displacement_vectors(&first[offset..], &second[offset..]);
        let b = displacement_vectors(&second[offset..], &third[offset..]);
        let c = displacement_vectors(&third[offset..], &fourth[offset..]);
        let left = cross(a, b);
        let right = cross(b, c);
        let axis_squared = dot(b, b).to_array();
        let left_squared = dot(left, left).to_array();
        let right_squared = dot(right, right).to_array();
        let signed = dot(cross(left, right), b).to_array();
        let cosine = dot(left, right).to_array();
        for lane in 0..LANES {
            output[offset + lane] = torsion_lane(
                axis_squared[lane],
                left_squared[lane],
                right_squared[lane],
                signed[lane],
                cosine[lane],
            );
        }
        offset += LANES;
    }
    for row in offset..first.len() {
        output[row] =
            match crate::measure::dihedral(first[row], second[row], third[row], fourth[row]) {
                Some(dihedral) => dihedral,
                None => f64::NAN,
            };
    }
}

fn torsion_lane(
    axis_squared: f64,
    left_squared: f64,
    right_squared: f64,
    signed: f64,
    cosine: f64,
) -> f64 {
    if axis_squared <= 0.0 || left_squared <= 0.0 || right_squared <= 0.0 {
        return f64::NAN;
    }
    let signed = signed * axis_squared.sqrt().recip();
    if matches!(signed.classify(), std::num::FpCategory::Zero)
        && matches!(cosine.classify(), std::num::FpCategory::Zero)
    {
        return f64::NAN;
    }
    signed.atan2(cosine)
}

fn displacement_vectors(from: &[[f32; 3]], to: &[[f32; 3]]) -> [F64xN; 3] {
    core::array::from_fn(|axis| {
        let from = F64xN::from(core::array::from_fn(|lane| f64::from(from[lane][axis])));
        let to = F64xN::from(core::array::from_fn(|lane| f64::from(to[lane][axis])));
        to - from
    })
}

fn dot(left: [F64xN; 3], right: [F64xN; 3]) -> F64xN {
    left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
}

fn cross(left: [F64xN; 3], right: [F64xN; 3]) -> [F64xN; 3] {
    [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    ]
}
