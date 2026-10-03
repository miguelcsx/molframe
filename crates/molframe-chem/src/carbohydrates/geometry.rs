//! Ring centroid, oriented polygon normal and anomeric direction.

use super::types::RingGeometry;
use molframe_core::{AtomIndex, Structure};
use num_traits::ToPrimitive;

pub(super) fn ring_geometry(
    structure: &Structure,
    ring: &[AtomIndex],
    anomeric: Option<AtomIndex>,
) -> Option<RingGeometry> {
    let mut positions = [[0.0_f64; 3]; 6];
    let mut center = [0.0; 3];
    for (index, atom) in ring.iter().enumerate() {
        let position = structure.atom(*atom)?.position()?;
        if !position.iter().all(|value| value.is_finite()) {
            return None;
        }
        for axis in 0..3 {
            positions[index][axis] = f64::from(position[axis]);
            center[axis] += positions[index][axis];
        }
    }
    let count = match ring.len() {
        5 => 5.0,
        6 => 6.0,
        _ => return None,
    };
    for value in &mut center {
        *value /= count;
    }
    let mut normal = [0.0; 3];
    for index in 0..ring.len() {
        let a = subtract(positions[index], center);
        let b = subtract(positions[(index + 1) % ring.len()], center);
        normal[0] += a[1] * b[2] - a[2] * b[1];
        normal[1] += a[2] * b[0] - a[0] * b[2];
        normal[2] += a[0] * b[1] - a[1] * b[0];
    }
    let anomeric_direction = match anomeric {
        Some(atom) => {
            let index = ring.iter().position(|value| *value == atom)?;
            unit(subtract(positions[index], center))
        }
        None => None,
    };
    Some(RingGeometry {
        center: [
            center[0].to_f32()?,
            center[1].to_f32()?,
            center[2].to_f32()?,
        ],
        normal: unit(normal)?,
        anomeric_direction,
    })
}

fn subtract(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn unit(vector: [f64; 3]) -> Option<[f32; 3]> {
    let length = vector.iter().map(|value| value * value).sum::<f64>().sqrt();
    if !length.is_finite() || length <= 1e-10 {
        return None;
    }
    Some([
        (vector[0] / length).to_f32()?,
        (vector[1] / length).to_f32()?,
        (vector[2] / length).to_f32()?,
    ])
}
