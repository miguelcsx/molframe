//! CHARMM CARTCV-style placement from three reference atoms.

use num_traits::ToPrimitive;
use pdbiox_geom::{cross, normalise};

/// Places `L` from references `I`, `J`, `K` and `K–L` geometry.
///
/// Returns `None` for a collinear reference frame or invalid scalar inputs.
#[must_use]
pub fn place_atom(
    i: [f32; 3],
    j: [f32; 3],
    k: [f32; 3],
    length: f64,
    angle: f64,
    torsion: f64,
) -> Option<[f32; 3]> {
    if !length.is_finite() || length < 0.0 || !angle.is_finite() || !torsion.is_finite() {
        return None;
    }
    let along = normalise(subtract(j, k))?;
    let incoming = subtract(i, j);
    let normal = normalise(cross(incoming, along))?;
    let plane = cross(along, normal);
    let radial = length * angle.sin();
    let vector = [
        length * angle.cos() * along[0]
            + radial * torsion.cos() * plane[0]
            + radial * torsion.sin() * normal[0],
        length * angle.cos() * along[1]
            + radial * torsion.cos() * plane[1]
            + radial * torsion.sin() * normal[1],
        length * angle.cos() * along[2]
            + radial * torsion.cos() * plane[2]
            + radial * torsion.sin() * normal[2],
    ];
    Some([
        (f64::from(k[0]) + vector[0]).to_f32()?,
        (f64::from(k[1]) + vector[1]).to_f32()?,
        (f64::from(k[2]) + vector[2]).to_f32()?,
    ])
}

fn subtract(left: [f32; 3], right: [f32; 3]) -> [f64; 3] {
    [
        f64::from(left[0] - right[0]),
        f64::from(left[1] - right[1]),
        f64::from(left[2] - right[2]),
    ]
}

#[cfg(test)]
#[path = "place_tests.rs"]
mod tests;
