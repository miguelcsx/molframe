//! Shared lattice-vector conversion for trajectory formats.

use pdbiox_core::structure::UnitCell;

pub(crate) fn cell_from_vectors(vectors: [[f64; 3]; 3]) -> Option<UnitCell> {
    let lengths = vectors.map(norm);
    if lengths.iter().any(|length| *length <= f64::EPSILON) {
        return None;
    }
    let angle = |left: [f64; 3], right: [f64; 3], left_length: f64, right_length: f64| {
        let cosine = dot(left, right) / (left_length * right_length);
        cosine.clamp(-1.0, 1.0).acos().to_degrees()
    };
    Some(UnitCell {
        lengths,
        angles: [
            angle(vectors[1], vectors[2], lengths[1], lengths[2]),
            angle(vectors[0], vectors[2], lengths[0], lengths[2]),
            angle(vectors[0], vectors[1], lengths[0], lengths[1]),
        ],
    })
}

pub(crate) fn vectors_from_cell(cell: UnitCell) -> Option<[[f64; 3]; 3]> {
    let [a, b, c] = cell.lengths;
    let [alpha, beta, gamma] = cell.angles.map(f64::to_radians);
    if [a, b, c, alpha, beta, gamma]
        .iter()
        .any(|value| !value.is_finite())
        || [a, b, c].iter().any(|length| *length <= 0.0)
    {
        return None;
    }
    let sin_gamma = gamma.sin();
    if sin_gamma.abs() <= f64::EPSILON {
        return None;
    }
    let cx = c * beta.cos();
    let cy = c * (alpha.cos() - beta.cos() * gamma.cos()) / sin_gamma;
    let cz_squared = c.mul_add(c, -cx.mul_add(cx, cy * cy));
    if cz_squared < -f64::EPSILON {
        return None;
    }
    Some([
        [a, 0.0, 0.0],
        [b * gamma.cos(), b * sin_gamma, 0.0],
        [cx, cy, cz_squared.max(0.0).sqrt()],
    ])
}

pub(crate) fn determinant(vectors: [[f64; 3]; 3]) -> f64 {
    vectors[0][0] * (vectors[1][1] * vectors[2][2] - vectors[1][2] * vectors[2][1])
        - vectors[0][1] * (vectors[1][0] * vectors[2][2] - vectors[1][2] * vectors[2][0])
        + vectors[0][2] * (vectors[1][0] * vectors[2][1] - vectors[1][1] * vectors[2][0])
}

fn dot(left: [f64; 3], right: [f64; 3]) -> f64 {
    left.iter().zip(right).map(|(a, b)| a * b).sum()
}

fn norm(vector: [f64; 3]) -> f64 {
    dot(vector, vector).sqrt()
}
