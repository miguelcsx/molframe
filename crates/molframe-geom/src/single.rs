//! Three-vector arithmetic in the precision positions are stored in (`f32`), for callers that
//! build a frame or a view from coordinates and have no use for a round trip through `f64`.
//!
//! The double-precision kernels in [`crate::measure`] are the ones to measure with. These are
//! plain component arithmetic, kept identical to what a caller would write inline.

/// `a - b`.
#[must_use]
#[inline]
pub fn sub_f32(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// The dot product of two vectors.
#[must_use]
#[inline]
pub fn dot_f32(a: [f32; 3], b: [f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// The cross product of two vectors.
#[must_use]
#[inline]
pub fn cross_f32(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// The vector scaled to unit length, or `None` when its length is not above `min_length` or is
/// not finite: a vector that short has no direction worth reporting.
#[must_use]
#[inline]
pub fn normalise_f32(vector: [f32; 3], min_length: f32) -> Option<[f32; 3]> {
    let length = dot_f32(vector, vector).sqrt();
    (length > min_length && length.is_finite())
        .then(|| [vector[0] / length, vector[1] / length, vector[2] / length])
}

#[cfg(test)]
#[path = "single_tests.rs"]
mod tests;
