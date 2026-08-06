//! Distances, angles and torsions.
//!
//! The kernels every other geometric quantity is built from. Two properties are
//! deliberate throughout:
//!
//! Squared distances are the primitive, and the square root is taken only when a
//! caller actually wants a length. A neighbour search comparing against a cutoff
//! never needs one, and a square root in that inner loop is the single most
//! common way a distance calculation ends up an order of magnitude slower than
//! it should be.
//!
//! Accumulation is in double precision even though positions are stored in
//! single. Storage precision is set by what was measured — three decimal places
//! in ångström — while summation precision is set by how many terms are summed,
//! and a centroid over a million atoms needs the wider type.

/// The vector from `from` to `to`.
#[must_use]
pub fn displacement(from: [f32; 3], to: [f32; 3]) -> [f64; 3] {
    [
        f64::from(to[0]) - f64::from(from[0]),
        f64::from(to[1]) - f64::from(from[1]),
        f64::from(to[2]) - f64::from(from[2]),
    ]
}

/// The squared distance between two positions.
///
/// Prefer this wherever the result is only compared: comparing squared
/// distances against a squared cutoff answers the same question without a
/// square root.
#[must_use]
pub fn distance_squared(a: [f32; 3], b: [f32; 3]) -> f64 {
    let d = displacement(a, b);
    d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
}

/// The distance between two positions.
///
/// # Examples
///
/// ```
/// use pdbiox_geom::distance;
///
/// let separation = distance([0.0, 0.0, 0.0], [3.0, 4.0, 0.0]);
/// assert!((separation - 5.0).abs() < 1e-12);
/// ```
#[must_use]
pub fn distance(a: [f32; 3], b: [f32; 3]) -> f64 {
    distance_squared(a, b).sqrt()
}

/// The dot product of two vectors.
#[must_use]
pub fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// The cross product of two vectors.
#[must_use]
pub fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// The length of a vector.
#[must_use]
pub fn norm(vector: [f64; 3]) -> f64 {
    dot(vector, vector).sqrt()
}

/// The vector scaled to unit length, or `None` when it has no direction.
///
/// A zero-length vector has no direction to report, and returning one anyway is
/// how a degenerate frame turns into a plausible-looking wrong angle.
#[must_use]
pub fn normalise(vector: [f64; 3]) -> Option<[f64; 3]> {
    let length = norm(vector);
    if length <= f64::EPSILON {
        return None;
    }
    Some([vector[0] / length, vector[1] / length, vector[2] / length])
}

/// The angle at `vertex` between `a` and `c`, in radians.
///
/// Returns `None` when either arm has no length, because the angle is then
/// undefined rather than zero.
///
/// # Examples
///
/// ```
/// use pdbiox_geom::angle;
///
/// let right = angle([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
/// assert!(right.is_some_and(|value| (value - std::f64::consts::FRAC_PI_2).abs() < 1e-12));
/// ```
#[must_use]
pub fn angle(a: [f32; 3], vertex: [f32; 3], c: [f32; 3]) -> Option<f64> {
    let first = normalise(displacement(vertex, a))?;
    let second = normalise(displacement(vertex, c))?;
    // Clamped because rounding can push the cosine a hair outside its range,
    // and the inverse cosine of 1.0000000001 is not a number.
    Some(dot(first, second).clamp(-1.0, 1.0).acos())
}

/// The torsion about the `b`–`c` bond, in radians, on the interval `(-π, π]`.
///
/// Returns `None` when the four positions do not define two planes — three
/// collinear points among them leave the torsion undefined.
///
/// # Examples
///
/// ```
/// use pdbiox_geom::dihedral;
///
/// let planar = dihedral(
///     [1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
/// );
/// assert!(planar.is_some_and(|value| value.abs() < 1e-9), "cis is zero");
/// ```
#[must_use]
pub fn dihedral(a: [f32; 3], b: [f32; 3], c: [f32; 3], d: [f32; 3]) -> Option<f64> {
    let first = displacement(a, b);
    let second = displacement(b, c);
    let third = displacement(c, d);

    let left = cross(first, second);
    let right = cross(second, third);
    let axis = normalise(second)?;

    // The signed angle between the two plane normals, taken about the central
    // bond so that the sign means the same thing on every torsion.
    let y = dot(cross(left, right), axis);
    let x = dot(left, right);
    if y == 0.0 && x == 0.0 {
        return None;
    }
    Some(y.atan2(x))
}

/// Converts radians to degrees, which is what torsions are reported in.
#[must_use]
pub fn degrees(radians: f64) -> f64 {
    radians.to_degrees()
}

#[cfg(test)]
#[path = "measure_tests.rs"]
mod tests;
