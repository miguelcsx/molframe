//! Deviation between two sets already in correspondence.
//!
//! These measure positions where they are; nothing is fitted. Use
//! [`superpose`](super::superpose) first if the sets need aligning.

use super::fit::{SuperposeError, exact_count};

/// The deviation between two sets already in correspondence.
///
/// No fitting happens: this measures the sets where they are. Use
/// [`superpose`](super::superpose) first if the sets need aligning.
///
/// Runs in `O(n)` time and `O(1)` auxiliary space.
///
/// # Errors
///
/// Returns [`SuperposeError::LengthMismatch`] when the sets differ in size.
pub fn rmsd(mobile: &[[f32; 3]], reference: &[[f32; 3]]) -> Result<f64, SuperposeError> {
    if mobile.len() != reference.len() {
        return Err(SuperposeError::LengthMismatch);
    }

    if mobile.is_empty() {
        return Ok(0.0);
    }

    let total = crate::simd::squared_deviation_sum(mobile, reference);

    let count = exact_count(mobile.len()).ok_or(SuperposeError::TooManyPoints)?;
    Ok((total / count).sqrt())
}

/// The deviation between two flat coordinate arrays.
///
/// The flat form exists for zero-copy callers holding one contiguous lane; the
/// values are otherwise identical to [`rmsd`].
///
/// Runs in `O(n)` time and `O(1)` auxiliary space.
///
/// # Errors
///
/// Returns [`SuperposeError::LengthMismatch`] when the arrays differ in length
/// or the length is not a whole number of three-component positions.
pub fn rmsd_flat(mobile: &[f32], reference: &[f32]) -> Result<f64, SuperposeError> {
    if mobile.len() != reference.len() || !mobile.len().is_multiple_of(3) {
        return Err(SuperposeError::LengthMismatch);
    }
    if mobile.is_empty() {
        return Ok(0.0);
    }

    let (mobile_blocks, _) = mobile.as_chunks::<3>();
    let (reference_blocks, _) = reference.as_chunks::<3>();
    let mut total = 0.0;
    for (left, right) in mobile_blocks.iter().zip(reference_blocks) {
        let dx = f64::from(left[0]) - f64::from(right[0]);
        let dy = f64::from(left[1]) - f64::from(right[1]);
        let dz = f64::from(left[2]) - f64::from(right[2]);
        total += dx.mul_add(dx, dy.mul_add(dy, dz * dz));
    }

    let count = exact_count(mobile.len() / 3).ok_or(SuperposeError::TooManyPoints)?;
    Ok((total / count).sqrt())
}

#[cfg(test)]
#[path = "deviation_tests.rs"]
mod tests;
