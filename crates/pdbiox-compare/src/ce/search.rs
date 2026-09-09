//! Allocation-reusing CE path search over fragment similarities.

use core::cmp::Ordering;

use super::workspace::{
    SearchWorkspace, bytes_for, checked_sum, empty_with_capacity, ensure_limit, memory_error,
    memory_plan,
};
use super::{CeError, CeOptions};
use crate::numeric::usize_to_f64;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct Fragment {
    pub(super) reference: usize,
    pub(super) mobile: usize,
}

#[derive(Debug, PartialEq)]
pub(super) struct Candidate {
    pub(super) fragments: Vec<Fragment>,
    pub(super) similarity: f64,
}

pub(super) fn find_paths(
    reference: &[[f32; 3]],
    mobile: &[[f32; 3]],
    options: CeOptions,
) -> Result<Vec<Candidate>, CeError> {
    let plan = memory_plan(reference.len(), mobile.len(), options)?;
    let workspace = SearchWorkspace::new(reference, mobile, options, &plan)?;
    let mut scratch = empty_with_capacity(plan.maximum_fragments, plan.minimum_bytes, options)?;
    let mut best = empty_with_capacity(options.max_paths, plan.minimum_bytes, options)?;
    let fixed_bytes =
        fixed_resident_bytes(&workspace, scratch.capacity(), best.capacity(), options)?;
    ensure_limit(fixed_bytes, options)?;
    let reference_starts = reference.len() - options.window_size + 1;
    let mobile_starts = mobile.len() - options.window_size + 1;
    for reference_start in 0..reference_starts {
        if cannot_beat(
            reference_start,
            reference.len(),
            &best,
            options.window_size,
            options.max_paths,
        ) {
            break;
        }
        for mobile_start in 0..mobile_starts {
            if cannot_beat(
                mobile_start,
                mobile.len(),
                &best,
                options.window_size,
                options.max_paths,
            ) {
                break;
            }
            let initial = workspace.get(reference_start, mobile_start);
            if initial <= options.fragment_similarity_threshold {
                continue;
            }
            scratch.clear();
            scratch.push(Fragment {
                reference: reference_start,
                mobile: mobile_start,
            });
            let similarity = extend_path(
                &mut scratch,
                initial,
                &workspace,
                reference,
                mobile,
                options,
            );
            retain_candidate(&mut best, &scratch, similarity, fixed_bytes, options)?;
        }
    }
    Ok(best)
}

fn fixed_resident_bytes(
    workspace: &SearchWorkspace,
    scratch_capacity: usize,
    best_capacity: usize,
    options: CeOptions,
) -> Result<usize, CeError> {
    let workspace_bytes = workspace.resident_bytes(options)?;
    let scratch_bytes = bytes_for::<Fragment>(scratch_capacity, options)?;
    let best_bytes = bytes_for::<Candidate>(best_capacity, options)?;
    checked_sum(
        checked_sum(workspace_bytes, scratch_bytes, options)?,
        best_bytes,
        options,
    )
}

pub(super) fn cannot_beat(
    start: usize,
    length: usize,
    best: &[Candidate],
    window: usize,
    max_paths: usize,
) -> bool {
    if best.len() < max_paths {
        return false;
    }
    best.last().is_some_and(|worst| {
        let gaps = worst.fragments.len().saturating_sub(1);
        let Some(span) = window.checked_mul(gaps) else {
            return true;
        };
        let Some(last_start) = length.checked_sub(span) else {
            return true;
        };
        start > last_start
    })
}

fn extend_path(
    fragments: &mut Vec<Fragment>,
    mut similarity: f64,
    workspace: &SearchWorkspace,
    reference: &[[f32; 3]],
    mobile: &[[f32; 3]],
    options: CeOptions,
) -> f64 {
    loop {
        if fragments.len() == fragments.capacity() {
            break;
        }
        let Some(last) = fragments.last().copied() else {
            break;
        };
        let Some(reference_base) = last.reference.checked_add(options.window_size) else {
            break;
        };
        let Some(mobile_base) = last.mobile.checked_add(options.window_size) else {
            break;
        };
        let reference_last = reference.len() - options.window_size;
        let mobile_last = mobile.len() - options.window_size;
        if reference_base > reference_last || mobile_base > mobile_last {
            break;
        }
        let useful_reference = reference_last - reference_base;
        let useful_mobile = mobile_last - mobile_base;
        let useful_gap_index = useful_reference
            .saturating_mul(2)
            .saturating_sub(1)
            .max(useful_mobile.saturating_mul(2));
        let maximum_gap_index = (options.max_gap * 2).min(useful_gap_index);
        let mut best_extension: Option<(Fragment, f64)> = None;
        for gap_index in 0..=maximum_gap_index {
            let mut next = Fragment {
                reference: reference_base,
                mobile: mobile_base,
            };
            if gap_index % 2 == 1 {
                next.reference += gap_index / 2 + 1;
            } else {
                next.mobile += gap_index / 2;
            }
            if next.reference > reference_last
                || next.mobile > mobile_last
                || workspace.get(next.reference, next.mobile)
                    <= options.fragment_similarity_threshold
            {
                continue;
            }
            let Some(cross) = path_cross_similarity(
                fragments,
                next,
                reference,
                mobile,
                options.window_size,
                options.path_similarity_threshold,
            ) else {
                continue;
            };
            if best_extension.is_none_or(|(_, best_cross)| cross > best_cross) {
                best_extension = Some((next, cross));
            }
        }
        let Some((next, cross)) = best_extension else {
            break;
        };
        let count = usize_to_f64(fragments.len());
        let current_terms = count + count * (count - 1.0) / 2.0;
        let new_terms = count + 1.0 + count * (count + 1.0) / 2.0;
        let next_similarity = (current_terms * similarity
            + count * cross
            + workspace.get(next.reference, next.mobile))
            / new_terms;
        if next_similarity <= options.path_similarity_threshold {
            break;
        }
        fragments.push(next);
        similarity = next_similarity;
    }
    similarity
}

fn path_cross_similarity(
    fragments: &[Fragment],
    next: Fragment,
    reference: &[[f32; 3]],
    mobile: &[[f32; 3]],
    window: usize,
    threshold: f64,
) -> Option<f64> {
    let count = usize_to_f64(fragments.len());
    let mut sum = 0.0;
    for &fragment in fragments {
        sum += cross_similarity(fragment, next, reference, mobile, window);
        let upper_bound = sum / count;
        if upper_bound <= threshold {
            return None;
        }
    }
    Some(sum / count)
}

fn cross_similarity(
    first: Fragment,
    second: Fragment,
    reference: &[[f32; 3]],
    mobile: &[[f32; 3]],
    window: usize,
) -> f64 {
    let mut difference =
        (pdbiox_geom::distance(reference[first.reference], reference[second.reference])
            - pdbiox_geom::distance(mobile[first.mobile], mobile[second.mobile]))
        .abs();
    difference += (pdbiox_geom::distance(
        reference[first.reference + window - 1],
        reference[second.reference + window - 1],
    ) - pdbiox_geom::distance(
        mobile[first.mobile + window - 1],
        mobile[second.mobile + window - 1],
    ))
    .abs();
    for offset in 1..window - 1 {
        difference += (pdbiox_geom::distance(
            reference[first.reference + offset],
            reference[second.reference + window - 1 - offset],
        ) - pdbiox_geom::distance(
            mobile[first.mobile + offset],
            mobile[second.mobile + window - 1 - offset],
        ))
        .abs();
    }
    -difference / usize_to_f64(window)
}

fn retain_candidate(
    best: &mut Vec<Candidate>,
    fragments: &[Fragment],
    similarity: f64,
    fixed_bytes: usize,
    options: CeOptions,
) -> Result<(), CeError> {
    let position = match best
        .iter()
        .position(|candidate| candidate_order(fragments.len(), similarity, candidate).is_lt())
    {
        Some(position) => position,
        None => best.len(),
    };
    if best.len() == options.max_paths && position == best.len() {
        return Ok(());
    }
    let mut storage = if best.len() == options.max_paths {
        match best.pop() {
            Some(candidate) => candidate.fragments,
            None => Vec::new(),
        }
    } else {
        Vec::new()
    };
    storage.clear();
    let other_bytes = retained_fragment_bytes(best, options)?;
    let requested_capacity = storage.capacity().max(fragments.len());
    let required = checked_sum(
        fixed_bytes,
        checked_sum(
            other_bytes,
            bytes_for::<Fragment>(requested_capacity, options)?,
            options,
        )?,
        options,
    )?;
    ensure_limit(required, options)?;
    if storage.capacity() < fragments.len() {
        storage
            .try_reserve_exact(fragments.len())
            .map_err(|_| memory_error(required, options))?;
    }
    let actual = checked_sum(
        fixed_bytes,
        checked_sum(
            other_bytes,
            bytes_for::<Fragment>(storage.capacity(), options)?,
            options,
        )?,
        options,
    )?;
    ensure_limit(actual, options)?;
    storage.extend_from_slice(fragments);
    best.insert(
        position,
        Candidate {
            fragments: storage,
            similarity,
        },
    );
    Ok(())
}

fn candidate_order(length: usize, similarity: f64, other: &Candidate) -> Ordering {
    other
        .fragments
        .len()
        .cmp(&length)
        .then_with(|| other.similarity.total_cmp(&similarity))
}

fn retained_fragment_bytes(candidates: &[Candidate], options: CeOptions) -> Result<usize, CeError> {
    let mut bytes = 0usize;
    for candidate in candidates {
        bytes = checked_sum(
            bytes,
            bytes_for::<Fragment>(candidate.fragments.capacity(), options)?,
            options,
        )?;
    }
    Ok(bytes)
}
