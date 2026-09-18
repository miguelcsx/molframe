//! Banded fragment geometry and the optional rectangular similarity cache.

use core::mem::size_of;

use super::search::{Candidate, Fragment};
use super::{CeError, CeOptions};
use crate::numeric::usize_to_f64;

pub(super) struct SearchMemoryPlan {
    pub(super) minimum_bytes: usize,
    pub(super) cache_bytes: Option<usize>,
    pub(super) maximum_fragments: usize,
    fixed_bytes: usize,
    worst_candidate_bytes: Option<usize>,
}

pub(super) fn memory_plan(
    reference_count: usize,
    mobile_count: usize,
    options: CeOptions,
) -> Result<SearchMemoryPlan, CeError> {
    let width = options.window_size - 2;
    let band_elements = checked_product(
        checked_sum(reference_count, mobile_count, options)?,
        width,
        options,
    )?;
    let band_bytes = bytes_for::<f64>(band_elements, options)?;
    let maximum_fragments = reference_count.min(mobile_count) / options.window_size;
    let scratch_bytes = bytes_for::<Fragment>(maximum_fragments, options)?;
    let best_bytes = bytes_for::<Candidate>(options.max_paths, options)?;
    let fixed_bytes = checked_sum(scratch_bytes, best_bytes, options)?;
    let minimum_bytes = checked_sum(band_bytes, fixed_bytes, options)?;
    ensure_limit(minimum_bytes, options)?;

    let worst_candidate_bytes = options
        .max_paths
        .checked_mul(maximum_fragments)
        .and_then(|count| count.checked_mul(size_of::<Fragment>()));
    let reference_starts = reference_count - options.window_size + 1;
    let mobile_starts = mobile_count - options.window_size + 1;
    let full_cache_bytes = reference_starts
        .checked_mul(mobile_starts)
        .and_then(|count| count.checked_mul(size_of::<f64>()));
    let cache_bytes = match (full_cache_bytes, worst_candidate_bytes) {
        (Some(cache), Some(candidates)) => minimum_bytes
            .checked_add(candidates)
            .and_then(|bytes| bytes.checked_add(cache))
            .filter(|bytes| *bytes <= options.memory_limit_bytes)
            .map(|_| cache),
        _ => None,
    };
    Ok(SearchMemoryPlan {
        minimum_bytes,
        cache_bytes,
        maximum_fragments,
        fixed_bytes,
        worst_candidate_bytes,
    })
}

pub(super) struct SearchWorkspace {
    reference: DistanceBand,
    mobile: DistanceBand,
    similarities: Option<SimilarityMatrix>,
    term_count: usize,
    fragment_similarity_threshold: f64,
}

impl SearchWorkspace {
    pub(super) fn new(
        reference: &[[f32; 3]],
        mobile: &[[f32; 3]],
        options: CeOptions,
        plan: &SearchMemoryPlan,
    ) -> Result<Self, CeError> {
        let reference_band =
            DistanceBand::new(reference, options.window_size, plan.minimum_bytes, options)?;
        let mobile_band =
            DistanceBand::new(mobile, options.window_size, plan.minimum_bytes, options)?;
        let term_count =
            checked_product(options.window_size - 1, options.window_size - 2, options)? / 2;
        let similarities = match plan.cache_bytes {
            Some(cache_bytes) => SimilarityMatrix::try_new(
                &reference_band,
                &mobile_band,
                term_count,
                cache_bytes,
                plan,
                options,
            ),
            None => None,
        };
        Ok(Self {
            reference: reference_band,
            mobile: mobile_band,
            similarities,
            term_count,
            fragment_similarity_threshold: options.fragment_similarity_threshold,
        })
    }

    pub(super) fn get(&self, reference_start: usize, mobile_start: usize) -> f64 {
        match &self.similarities {
            Some(similarities) => similarities.get(reference_start, mobile_start),
            None => fragment_similarity(
                &self.reference,
                &self.mobile,
                reference_start,
                mobile_start,
                self.term_count,
                self.fragment_similarity_threshold,
            ),
        }
    }

    pub(super) fn resident_bytes(&self, options: CeOptions) -> Result<usize, CeError> {
        let bands = checked_sum(
            bytes_for::<f64>(self.reference.values.capacity(), options)?,
            bytes_for::<f64>(self.mobile.values.capacity(), options)?,
            options,
        )?;
        match &self.similarities {
            Some(similarities) => checked_sum(
                bands,
                bytes_for::<f64>(similarities.values.capacity(), options)?,
                options,
            ),
            None => Ok(bands),
        }
    }
}

struct DistanceBand {
    width: usize,
    values: Vec<f64>,
}

impl DistanceBand {
    fn new(
        points: &[[f32; 3]],
        window: usize,
        required: usize,
        options: CeOptions,
    ) -> Result<Self, CeError> {
        let width = window - 2;
        let count = checked_product(points.len(), width, options)?;
        let mut values = empty_with_capacity(count, required, options)?;
        values.resize(count, 0.0);
        for first in 0..points.len() {
            let maximum_delta = (points.len() - first - 1).min(window - 1);
            for delta in 2..=maximum_delta {
                values[first * width + delta - 2] =
                    molframe_geom::distance(points[first], points[first + delta]);
            }
        }
        Ok(Self { width, values })
    }
}

struct SimilarityMatrix {
    columns: usize,
    values: Vec<f64>,
}

impl SimilarityMatrix {
    fn try_new(
        reference: &DistanceBand,
        mobile: &DistanceBand,
        term_count: usize,
        cache_bytes: usize,
        plan: &SearchMemoryPlan,
        options: CeOptions,
    ) -> Option<Self> {
        let columns = mobile.values.len() / mobile.width - options.window_size + 1;
        let rows = reference.values.len() / reference.width - options.window_size + 1;
        let count = cache_bytes / size_of::<f64>();
        let mut values = Vec::new();
        if values.try_reserve_exact(count).is_err() {
            return None;
        }
        let actual_cache = values.capacity().checked_mul(size_of::<f64>())?;
        let bands = reference
            .values
            .capacity()
            .checked_add(mobile.values.capacity())?
            .checked_mul(size_of::<f64>())?;
        let candidates = plan.worst_candidate_bytes?;
        let peak = bands
            .checked_add(plan.fixed_bytes)?
            .checked_add(candidates)?
            .checked_add(actual_cache)?;
        if peak > options.memory_limit_bytes {
            return None;
        }
        for reference_start in 0..rows {
            for mobile_start in 0..columns {
                values.push(fragment_similarity(
                    reference,
                    mobile,
                    reference_start,
                    mobile_start,
                    term_count,
                    options.fragment_similarity_threshold,
                ));
            }
        }
        Some(Self { columns, values })
    }

    fn get(&self, row: usize, column: usize) -> f64 {
        self.values[row * self.columns + column]
    }
}

fn fragment_similarity(
    reference: &DistanceBand,
    mobile: &DistanceBand,
    reference_start: usize,
    mobile_start: usize,
    term_count: usize,
    rejection_threshold: f64,
) -> f64 {
    let mut difference = 0.0;
    let window = reference.width + 2;
    let normalizer = usize_to_f64(term_count);
    for left in 0..window - 2 {
        let count = window - left - 2;
        let reference_offset = (reference_start + left) * reference.width;
        let mobile_offset = (mobile_start + left) * mobile.width;
        let reference_terms = &reference.values[reference_offset..reference_offset + count];
        let mobile_terms = &mobile.values[mobile_offset..mobile_offset + count];
        for (&reference_distance, &mobile_distance) in reference_terms.iter().zip(mobile_terms) {
            difference += (reference_distance - mobile_distance).abs();
        }
        if left + 1 < window - 2 && -difference / normalizer <= rejection_threshold {
            return rejection_threshold;
        }
    }
    -difference / normalizer
}

pub(super) fn empty_with_capacity<T>(
    capacity: usize,
    required: usize,
    options: CeOptions,
) -> Result<Vec<T>, CeError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(capacity)
        .map_err(|_| memory_error(required, options))?;
    Ok(values)
}

pub(super) fn bytes_for<T>(count: usize, options: CeOptions) -> Result<usize, CeError> {
    checked_product(count, size_of::<T>(), options)
}

pub(super) fn checked_product(
    first: usize,
    second: usize,
    options: CeOptions,
) -> Result<usize, CeError> {
    first
        .checked_mul(second)
        .ok_or(memory_error(usize::MAX, options))
}

pub(super) fn checked_sum(
    first: usize,
    second: usize,
    options: CeOptions,
) -> Result<usize, CeError> {
    first
        .checked_add(second)
        .ok_or(memory_error(usize::MAX, options))
}

pub(super) fn ensure_limit(required: usize, options: CeOptions) -> Result<(), CeError> {
    if required > options.memory_limit_bytes {
        Err(memory_error(required, options))
    } else {
        Ok(())
    }
}

pub(super) const fn memory_error(required: usize, options: CeOptions) -> CeError {
    CeError::MemoryLimit {
        required,
        limit: options.memory_limit_bytes,
    }
}
