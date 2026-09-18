//! Bounded candidate expansion and allocation-free final ordering.

use core::cmp::Ordering;

use super::search::{Candidate, Fragment};
use super::workspace::{
    bytes_for, checked_product, checked_sum, empty_with_capacity, ensure_limit,
};
use super::{CeAlignment, CeError, CeOptions, significance};

pub(super) fn finish_candidates(
    candidates: Vec<Candidate>,
    reference: &[[f32; 3]],
    mobile: &[[f32; 3]],
    options: CeOptions,
) -> Result<Vec<CeAlignment>, CeError> {
    let (required, maximum_points) = finish_memory(&candidates, candidates.capacity(), options)?;
    ensure_limit(required, options)?;
    let mut points = PointWorkspace {
        reference: empty_with_capacity(maximum_points, required, options)?,
        mobile: empty_with_capacity(maximum_points, required, options)?,
    };
    let mut alignments = empty_with_capacity(candidates.len(), required, options)?;
    for candidate in candidates {
        alignments.push(finish_one(
            &candidate,
            reference,
            mobile,
            options,
            required,
            &mut points,
        )?);
    }
    stable_insertion_sort(&mut alignments);
    Ok(alignments)
}

fn finish_memory(
    candidates: &[Candidate],
    candidate_capacity: usize,
    options: CeOptions,
) -> Result<(usize, usize), CeError> {
    let mut candidate_storage = bytes_for::<Candidate>(candidate_capacity, options)?;
    let mut output_indices = 0usize;
    let mut maximum_points = 0usize;
    for candidate in candidates {
        candidate_storage = checked_sum(
            candidate_storage,
            bytes_for::<Fragment>(candidate.fragments.capacity(), options)?,
            options,
        )?;
        let point_count = checked_product(candidate.fragments.len(), options.window_size, options)?;
        maximum_points = maximum_points.max(point_count);
        output_indices = checked_sum(
            output_indices,
            checked_product(bytes_for::<usize>(point_count, options)?, 2, options)?,
            options,
        )?;
    }
    let alignment_storage = bytes_for::<CeAlignment>(candidates.len(), options)?;
    let point_workspace =
        checked_product(bytes_for::<[f32; 3]>(maximum_points, options)?, 2, options)?;
    let required = checked_sum(
        checked_sum(candidate_storage, output_indices, options)?,
        checked_sum(alignment_storage, point_workspace, options)?,
        options,
    )?;
    Ok((required, maximum_points))
}

struct PointWorkspace {
    reference: Vec<[f32; 3]>,
    mobile: Vec<[f32; 3]>,
}

fn finish_one(
    candidate: &Candidate,
    reference: &[[f32; 3]],
    mobile: &[[f32; 3]],
    options: CeOptions,
    required: usize,
    points: &mut PointWorkspace,
) -> Result<CeAlignment, CeError> {
    let point_count = checked_product(candidate.fragments.len(), options.window_size, options)?;
    let mut reference_indices = empty_with_capacity(point_count, required, options)?;
    let mut mobile_indices = empty_with_capacity(point_count, required, options)?;
    points.reference.clear();
    points.mobile.clear();
    for fragment in &candidate.fragments {
        for offset in 0..options.window_size {
            let reference_index = fragment.reference + offset;
            let mobile_index = fragment.mobile + offset;
            reference_indices.push(reference_index);
            mobile_indices.push(mobile_index);
            points.reference.push(reference[reference_index]);
            points.mobile.push(mobile[mobile_index]);
        }
    }
    let rmsd = molframe_geom::superpose(&points.mobile, &points.reference)
        .map_err(CeError::Superpose)?
        .rmsd;
    let z_score = options
        .significance
        .map(|profile| significance(candidate, profile));
    Ok(CeAlignment {
        reference_indices,
        mobile_indices,
        fragment_count: candidate.fragments.len(),
        similarity: candidate.similarity,
        z_score,
        rmsd,
    })
}

fn stable_insertion_sort(alignments: &mut [CeAlignment]) {
    for index in 1..alignments.len() {
        let mut current = index;
        while current > 0 && alignment_order(&alignments[current], &alignments[current - 1]).is_lt()
        {
            alignments.swap(current, current - 1);
            current -= 1;
        }
    }
}

fn alignment_order(left: &CeAlignment, right: &CeAlignment) -> Ordering {
    right
        .reference_indices
        .len()
        .cmp(&left.reference_indices.len())
        .then_with(|| left.rmsd.total_cmp(&right.rmsd))
        .then_with(|| right.similarity.total_cmp(&left.similarity))
}
