//! Paired centres, scatter tensors and cross-covariance in one pass.
//!
//! Split from the fitting itself so each file stays within the size ceiling.
//!
//! The point-at-a-time Welford recurrence is numerically excellent and
//! structurally poor: every point needs a reciprocal, and each update reads the
//! centre the previous update wrote, so the loop is a chain of dependent
//! divisions that neither pipelines nor vectorises.
//!
//! Blocking recovers both. A block of fixed size is reduced with plain sums,
//! which have no loop-carried division, and blocks are combined with the
//! parallel-merge form of the same recurrence. Accuracy stays close to the
//! online form because no block sums more than `BLOCK_POINTS` deviations, and
//! the block size is a constant rather than a worker count, so the result is
//! identical on every machine (FR-515).

/// Points reduced per block before merging.
///
/// Small enough that a block's plain sums cannot drift, large enough that the
/// per-block merge and its six divisions disappear against the block itself.
const BLOCK_POINTS: usize = 64;

/// Statistics required to solve and validate a rigid fit.
#[derive(Clone, Copy)]
pub(super) struct FitStatistics {
    pub(super) mobile_centre: [f64; 3],
    pub(super) reference_centre: [f64; 3],
    pub(super) covariance: [[f64; 3]; 3],
    pub(super) mobile_scatter: [[f64; 3]; 3],
    pub(super) reference_scatter: [[f64; 3]; 3],
}

/// Computes paired centres, scatter tensors and cross-covariance in one pass.
///
/// Runs in `O(n)` time and `O(1)` space. Points beyond the shorter of the two
/// sets are ignored, which is what pairing them means.
pub(super) fn fit_statistics(mobile: &[[f32; 3]], reference: &[[f32; 3]]) -> FitStatistics {
    let paired = mobile.len().min(reference.len());
    let (Some(mobile), Some(reference)) = (mobile.get(..paired), reference.get(..paired)) else {
        return Moments::default().finish();
    };

    let mut total = Moments::default();
    for (mobile_block, reference_block) in mobile
        .chunks(BLOCK_POINTS)
        .zip(reference.chunks(BLOCK_POINTS))
    {
        total.merge(&Moments::of_block(mobile_block, reference_block));
    }

    total.finish()
}

/// Centres and second moments of some prefix of the paired sets.
#[derive(Clone, Copy, Default)]
struct Moments {
    /// How many pairs are summarised.
    count: f64,
    /// Mean of the mobile positions.
    mobile_centre: [f64; 3],
    /// Mean of the reference positions.
    reference_centre: [f64; 3],
    /// Co-moment of the two centred sets.
    covariance: [[f64; 3]; 3],
    /// Second moment of the centred mobile set.
    mobile_scatter: [[f64; 3]; 3],
    /// Second moment of the centred reference set.
    reference_scatter: [[f64; 3]; 3],
}

impl Moments {
    /// Reduces one block with plain sums and a single division per axis.
    fn of_block(mobile: &[[f32; 3]], reference: &[[f32; 3]]) -> Self {
        let count = exact_len(mobile.len().min(reference.len()));
        if count == 0.0 {
            return Self::default();
        }

        let mut mobile_centre = crate::simd::sum_positions(mobile);
        let mut reference_centre = crate::simd::sum_positions(reference);
        let inverse = count.recip();
        for axis in 0..3 {
            mobile_centre[axis] *= inverse;
            reference_centre[axis] *= inverse;
        }

        let mut covariance = [[0.0f64; 3]; 3];
        let mut mobile_scatter = [[0.0f64; 3]; 3];
        let mut reference_scatter = [[0.0f64; 3]; 3];
        crate::simd::paired_moments(
            mobile,
            reference,
            mobile_centre,
            reference_centre,
            &mut covariance,
            &mut mobile_scatter,
            &mut reference_scatter,
        );

        Self {
            count,
            mobile_centre,
            reference_centre,
            covariance,
            mobile_scatter,
            reference_scatter,
        }
    }

    /// Combines `other`, which covers the pairs immediately after these.
    ///
    /// This is the parallel-merge form of the Welford recurrence: the moments
    /// add, plus a correction for the two blocks' centres being apart.
    fn merge(&mut self, other: &Self) {
        let total = self.count + other.count;
        if total == 0.0 {
            return;
        }

        let mobile_gap = subtract(other.mobile_centre, self.mobile_centre);
        let reference_gap = subtract(other.reference_centre, self.reference_centre);
        let correction = self.count * other.count / total;
        let share = other.count / total;

        for row in 0..3 {
            for column in 0..3 {
                self.covariance[row][column] += other.covariance[row][column]
                    + mobile_gap[row] * reference_gap[column] * correction;
                self.mobile_scatter[row][column] += other.mobile_scatter[row][column]
                    + mobile_gap[row] * mobile_gap[column] * correction;
                self.reference_scatter[row][column] += other.reference_scatter[row][column]
                    + reference_gap[row] * reference_gap[column] * correction;
            }
        }
        for axis in 0..3 {
            self.mobile_centre[axis] += mobile_gap[axis] * share;
            self.reference_centre[axis] += reference_gap[axis] * share;
        }
        self.count = total;
    }

    /// Drops the running count, which the fit does not need.
    fn finish(self) -> FitStatistics {
        FitStatistics {
            mobile_centre: self.mobile_centre,
            reference_centre: self.reference_centre,
            covariance: self.covariance,
            mobile_scatter: self.mobile_scatter,
            reference_scatter: self.reference_scatter,
        }
    }
}

/// Computes weighted paired centres, scatter tensors and cross-covariance in one
/// pass.
///
/// An empty `weights` weights every pair one, which reproduces [`fit_statistics`]
/// exactly. A supplied slice is read position by position; a shorter slice leaves
/// the remaining pairs at weight one, and validation lives with the caller.
///
/// Runs in `O(n)` time and `O(1)` space. Points beyond the shorter of the two
/// sets are ignored, which is what pairing them means.
pub(super) fn fit_statistics_weighted(
    mobile: &[[f32; 3]],
    reference: &[[f32; 3]],
    weights: &[f64],
) -> FitStatistics {
    if weights.is_empty() {
        return fit_statistics(mobile, reference);
    }
    let paired_len = mobile.len().min(reference.len());
    let (Some(mobile), Some(reference)) = (mobile.get(..paired_len), reference.get(..paired_len))
    else {
        return Moments::default().finish();
    };
    // A short weight slice leaves the pairs past its end at weight one; the
    // fallback lives in `of_block` so no block ever silently loses its tail.
    let mut total = WeightedMoments::default();
    for start in (0..paired_len).step_by(BLOCK_POINTS) {
        let end = (start + BLOCK_POINTS).min(paired_len);
        let (Some(mobile_block), Some(reference_block)) =
            (mobile.get(start..end), reference.get(start..end))
        else {
            continue;
        };
        total.merge(&WeightedMoments::of_block(
            mobile_block,
            reference_block,
            &WeightBlock {
                weights,
                offset: start,
            },
        ));
    }
    total.finish()
}

/// A view of the caller's weight slice starting at one pair's position.
///
/// Reading through this rather than slicing means a weight list shorter than the
/// paired sets leaves the pairs past its end at weight one, instead of the block
/// silently shrinking to the shorter of the two — which would drop the tail of
/// every structure whose weight list is short.
#[derive(Clone, Copy)]
struct WeightBlock<'a> {
    weights: &'a [f64],
    offset: usize,
}

impl WeightBlock<'_> {
    fn of(&self, position: usize) -> f64 {
        let absolute = self.offset + position;
        match self.weights.get(absolute) {
            Some(weight) => *weight,
            None => 1.0,
        }
    }
}

/// Centres and weighted second moments of some prefix of the paired sets.
#[derive(Clone, Copy, Default)]
struct WeightedMoments {
    /// Total pair weight summarised.
    weight: f64,
    mobile_centre: [f64; 3],
    reference_centre: [f64; 3],
    covariance: [[f64; 3]; 3],
    mobile_scatter: [[f64; 3]; 3],
    reference_scatter: [[f64; 3]; 3],
}

impl WeightedMoments {
    fn of_block(mobile: &[[f32; 3]], reference: &[[f32; 3]], weights: &WeightBlock<'_>) -> Self {
        let mut weight = 0.0f64;
        let mut mobile_centre = [0.0f64; 3];
        let mut reference_centre = [0.0f64; 3];
        for (position, (&mobile_point, &reference_point)) in
            mobile.iter().zip(reference).enumerate()
        {
            let pair_weight = weights.of(position);
            weight += pair_weight;
            for axis in 0..3 {
                mobile_centre[axis] += pair_weight * f64::from(mobile_point[axis]);
                reference_centre[axis] += pair_weight * f64::from(reference_point[axis]);
            }
        }
        // A total weight that is not positive has no centre to divide by, and a
        // NaN would propagate through every term.
        if weight.is_nan() || weight <= 0.0 {
            return Self::default();
        }
        let inverse = weight.recip();
        for axis in 0..3 {
            mobile_centre[axis] *= inverse;
            reference_centre[axis] *= inverse;
        }

        let mut covariance = [[0.0f64; 3]; 3];
        let mut mobile_scatter = [[0.0f64; 3]; 3];
        let mut reference_scatter = [[0.0f64; 3]; 3];
        for (position, (&mobile_point, &reference_point)) in
            mobile.iter().zip(reference).enumerate()
        {
            let pair_weight = weights.of(position);
            let mobile_offset: [f64; 3] =
                std::array::from_fn(|axis| f64::from(mobile_point[axis]) - mobile_centre[axis]);
            let reference_offset: [f64; 3] = std::array::from_fn(|axis| {
                f64::from(reference_point[axis]) - reference_centre[axis]
            });
            for row in 0..3 {
                for column in 0..3 {
                    covariance[row][column] +=
                        pair_weight * mobile_offset[row] * reference_offset[column];
                    mobile_scatter[row][column] +=
                        pair_weight * mobile_offset[row] * mobile_offset[column];
                    reference_scatter[row][column] +=
                        pair_weight * reference_offset[row] * reference_offset[column];
                }
            }
        }

        Self {
            weight,
            mobile_centre,
            reference_centre,
            covariance,
            mobile_scatter,
            reference_scatter,
        }
    }

    /// Combines `other`, which covers the pairs immediately after these.
    fn merge(&mut self, other: &Self) {
        let total = self.weight + other.weight;
        if total <= 0.0 {
            return;
        }
        let mobile_gap = subtract(other.mobile_centre, self.mobile_centre);
        let reference_gap = subtract(other.reference_centre, self.reference_centre);
        let correction = self.weight * other.weight / total;
        let share = other.weight / total;

        for row in 0..3 {
            for column in 0..3 {
                self.covariance[row][column] += other.covariance[row][column]
                    + mobile_gap[row] * reference_gap[column] * correction;
                self.mobile_scatter[row][column] += other.mobile_scatter[row][column]
                    + mobile_gap[row] * mobile_gap[column] * correction;
                self.reference_scatter[row][column] += other.reference_scatter[row][column]
                    + reference_gap[row] * reference_gap[column] * correction;
            }
        }
        for axis in 0..3 {
            self.mobile_centre[axis] += mobile_gap[axis] * share;
            self.reference_centre[axis] += reference_gap[axis] * share;
        }
        self.weight = total;
    }

    fn finish(self) -> FitStatistics {
        FitStatistics {
            mobile_centre: self.mobile_centre,
            reference_centre: self.reference_centre,
            covariance: self.covariance,
            mobile_scatter: self.mobile_scatter,
            reference_scatter: self.reference_scatter,
        }
    }
}

/// Subtracts two three-dimensional vectors.
#[inline]
fn subtract(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [left[0] - right[0], left[1] - right[1], left[2] - right[2]]
}

/// A block length as a double, exactly for every block this code produces.
#[inline]
fn exact_len(count: usize) -> f64 {
    match u32::try_from(count) {
        Ok(count) => f64::from(count),
        Err(_) => f64::INFINITY,
    }
}

#[cfg(test)]
#[path = "fit_statistics_tests.rs"]
mod tests;
