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
///
/// Runs in `O(1)` time, performs all subtraction in `f64`, and allocates no
/// heap memory.
#[must_use]
#[inline]
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
///
/// Runs in `O(1)` time and allocates no memory.
#[must_use]
#[inline]
pub fn distance_squared(a: [f32; 3], b: [f32; 3]) -> f64 {
    let dx = f64::from(b[0]) - f64::from(a[0]);
    let dy = f64::from(b[1]) - f64::from(a[1]);
    let dz = f64::from(b[2]) - f64::from(a[2]);

    dx * dx + dy * dy + dz * dz
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
///
/// Runs in `O(1)` time and performs exactly one square root.
#[must_use]
#[inline]
pub fn distance(a: [f32; 3], b: [f32; 3]) -> f64 {
    distance_squared(a, b).sqrt()
}

/// The dot product of two vectors.
///
/// Runs in `O(1)` time and allocates no memory.
#[must_use]
#[inline]
pub fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// The cross product of two vectors.
///
/// Runs in `O(1)` time and allocates no heap memory.
#[must_use]
#[inline]
pub fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// The length of a vector.
///
/// Runs in `O(1)` time and performs exactly one square root.
#[must_use]
#[inline]
pub fn norm(vector: [f64; 3]) -> f64 {
    norm_squared(vector).sqrt()
}

/// The vector scaled to unit length, or `None` when it has no direction.
///
/// A zero-length vector has no direction to report, and returning one anyway is
/// how a degenerate frame turns into a plausible-looking wrong angle.
///
/// Runs in `O(1)` time, performs one square root and one division, and allocates
/// no heap memory.
#[must_use]
#[inline]
pub fn normalise(vector: [f64; 3]) -> Option<[f64; 3]> {
    let squared = norm_squared(vector);
    if squared <= 0.0 {
        return None;
    }

    let inverse_length = squared.sqrt().recip();
    Some([
        vector[0] * inverse_length,
        vector[1] * inverse_length,
        vector[2] * inverse_length,
    ])
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
///
/// Runs in `O(1)` time and allocates no heap memory.
#[must_use]
pub fn angle(a: [f32; 3], vertex: [f32; 3], c: [f32; 3]) -> Option<f64> {
    let first = displacement(vertex, a);
    let second = displacement(vertex, c);
    let first_squared = norm_squared(first);
    let second_squared = norm_squared(second);

    if first_squared <= 0.0 || second_squared <= 0.0 {
        return None;
    }

    let inverse_lengths = (first_squared * second_squared).sqrt().recip();

    // Clamped because rounding can push the cosine a hair outside its range,
    // and the inverse cosine of 1.0000000001 is not a number.
    Some(
        (dot(first, second) * inverse_lengths)
            .clamp(-1.0, 1.0)
            .acos(),
    )
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
///
/// Runs in `O(1)` time and allocates no heap memory.
#[must_use]
pub fn dihedral(
    first_point: [f32; 3],
    second_point: [f32; 3],
    third_point: [f32; 3],
    fourth_point: [f32; 3],
) -> Option<f64> {
    let first = displacement(first_point, second_point);
    let second = displacement(second_point, third_point);
    let third = displacement(third_point, fourth_point);

    let left = cross(first, second);
    let right = cross(second, third);
    let second_squared = norm_squared(second);

    if second_squared <= 0.0 || norm_squared(left) <= 0.0 || norm_squared(right) <= 0.0 {
        return None;
    }

    // The signed angle between the two plane normals, taken about the central
    // bond so that the sign means the same thing on every torsion.
    let inverse_axis_length = second_squared.sqrt().recip();
    let signed_component = dot(cross(left, right), second) * inverse_axis_length;
    let cosine_component = dot(left, right);

    if matches!(signed_component.classify(), std::num::FpCategory::Zero)
        && matches!(cosine_component.classify(), std::num::FpCategory::Zero)
    {
        return None;
    }

    Some(signed_component.atan2(cosine_component))
}

/// Converts radians to degrees, which is what torsions are reported in.
///
/// Runs in `O(1)` time and allocates no memory.
#[must_use]
#[inline]
pub fn degrees(radians: f64) -> f64 {
    radians.to_degrees()
}

/// Returns the squared length of a vector.
///
/// This is the internal primitive used when only a zero test or a later single
/// square root is required. Runs in `O(1)` time and allocates no memory.
#[inline]
fn norm_squared(vector: [f64; 3]) -> f64 {
    dot(vector, vector)
}

#[cfg(test)]
#[path = "measure_tests.rs"]
mod tests;
