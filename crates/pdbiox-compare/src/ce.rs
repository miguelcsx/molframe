//! Combinatorial Extension (CE) structural alignment over guide coordinates.

use pdbiox_geom::SuperposeError;

use crate::numeric::usize_to_f64;

/// Named significance calibration for a CE candidate path.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CeSignificanceProfile {
    /// Original CE calibration for eight-residue fragments.
    OriginalWindowEight,
}

/// Search controls for CE aligned-fragment paths.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CeOptions {
    /// Number of consecutive guide atoms in one aligned fragment pair.
    pub window_size: usize,
    /// Maximum insertion examined between consecutive fragments.
    pub max_gap: usize,
    /// Number of highest-ranked candidate paths retained.
    pub max_paths: usize,
    /// Minimum fragment-pair similarity; values nearer zero are better.
    pub fragment_similarity_threshold: f64,
    /// Minimum cross-fragment and accumulated-path similarity.
    pub path_similarity_threshold: f64,
    /// Optional named statistical calibration for reported significance.
    pub significance: Option<CeSignificanceProfile>,
}

impl CeOptions {
    /// Original CE search and significance profile.
    #[must_use]
    pub const fn original() -> Self {
        Self {
            window_size: 8,
            max_gap: 30,
            max_paths: 20,
            fragment_similarity_threshold: -3.0,
            path_similarity_threshold: -4.0,
            significance: Some(CeSignificanceProfile::OriginalWindowEight),
        }
    }
}

/// One CE correspondence with its structural fit.
#[derive(Clone, Debug, PartialEq)]
pub struct CeAlignment {
    /// Indices into the reference guide coordinates.
    pub reference_indices: Vec<usize>,
    /// Corresponding indices into the mobile guide coordinates.
    pub mobile_indices: Vec<usize>,
    /// Number of aligned fragment pairs before expansion to guide atoms.
    pub fragment_count: usize,
    /// CE path similarity, where values nearer zero are better.
    pub similarity: f64,
    /// Empirical CE significance estimate when a named profile was selected.
    pub z_score: Option<f64>,
    /// RMSD after shared rigid superposition of the correspondence.
    pub rmsd: f64,
}

/// Invalid CE settings, inputs or final fit.
#[derive(Clone, Copy, Debug, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum CeError {
    /// Window, gap or path count is outside the supported domain.
    #[error("invalid CE search options")]
    InvalidOptions,
    /// A coordinate set cannot contain two complete fragments.
    #[error("CE needs at least {required} guide atoms, found {actual}")]
    TooFewPoints {
        /// Minimum coordinate count for these options.
        required: usize,
        /// Supplied coordinate count.
        actual: usize,
    },
    /// No fragment pair passed the CE similarity threshold.
    #[error("CE found no compatible aligned fragment path")]
    NoAlignment,
    /// Shared rigid superposition failed.
    #[error("CE superposition failed: {0:?}")]
    Superpose(SuperposeError),
}

/// Finds the best CE alignment, preferring maximum coverage then minimum RMSD.
///
/// # Errors
///
/// Returns an option, input, no-path or superposition error.
pub fn ce_align(
    reference: &[[f32; 3]],
    mobile: &[[f32; 3]],
    options: CeOptions,
) -> Result<CeAlignment, CeError> {
    ce_alignments(reference, mobile, options)?
        .into_iter()
        .next()
        .ok_or(CeError::NoAlignment)
}

/// Returns the retained CE candidates ordered by coverage and fitted RMSD.
///
/// Runtime is `O(R² + M² + R M P (G + P W))` in guide atoms, retained path
/// length `P`, gap bound `G` and window `W`; memory is `O(R² + M² + R M)`.
///
/// # Errors
///
/// Returns an option, input, no-path or superposition error.
pub fn ce_alignments(
    reference: &[[f32; 3]],
    mobile: &[[f32; 3]],
    options: CeOptions,
) -> Result<Vec<CeAlignment>, CeError> {
    validate(reference, mobile, options)?;
    let reference_distances = DistanceMatrix::new(reference);
    let mobile_distances = DistanceMatrix::new(mobile);
    let similarity =
        similarity_matrix(&reference_distances, &mobile_distances, options.window_size);
    let candidates = find_paths(
        &similarity,
        &reference_distances,
        &mobile_distances,
        options,
    );
    if candidates.is_empty() {
        return Err(CeError::NoAlignment);
    }
    let mut alignments = candidates
        .into_iter()
        .map(|candidate| finish(&candidate, reference, mobile, options))
        .collect::<Result<Vec<_>, _>>()?;
    alignments.sort_by(|left, right| {
        right
            .reference_indices
            .len()
            .cmp(&left.reference_indices.len())
            .then_with(|| left.rmsd.total_cmp(&right.rmsd))
            .then_with(|| right.similarity.total_cmp(&left.similarity))
    });
    Ok(alignments)
}

fn validate(
    reference: &[[f32; 3]],
    mobile: &[[f32; 3]],
    options: CeOptions,
) -> Result<(), CeError> {
    if options.window_size < 3
        || options.max_paths == 0
        || options.max_gap.checked_mul(2).is_none()
        || !options.fragment_similarity_threshold.is_finite()
        || !options.path_similarity_threshold.is_finite()
        || options.fragment_similarity_threshold <= options.path_similarity_threshold
        || (matches!(
            options.significance,
            Some(CeSignificanceProfile::OriginalWindowEight)
        ) && options.window_size != 8)
    {
        return Err(CeError::InvalidOptions);
    }
    let required = options
        .window_size
        .checked_mul(2)
        .ok_or(CeError::InvalidOptions)?;
    for actual in [reference.len(), mobile.len()] {
        if actual < required {
            return Err(CeError::TooFewPoints { required, actual });
        }
    }
    if reference
        .iter()
        .chain(mobile)
        .flatten()
        .any(|value| !value.is_finite())
    {
        return Err(CeError::InvalidOptions);
    }
    Ok(())
}

struct DistanceMatrix {
    size: usize,
    values: Vec<f64>,
}

impl DistanceMatrix {
    fn new(points: &[[f32; 3]]) -> Self {
        let mut values = vec![0.0; points.len() * points.len()];
        for row in 0..points.len() {
            for column in row + 1..points.len() {
                let distance = pdbiox_geom::distance(points[row], points[column]);
                values[row * points.len() + column] = distance;
                values[column * points.len() + row] = distance;
            }
        }
        Self {
            size: points.len(),
            values,
        }
    }

    fn get(&self, row: usize, column: usize) -> f64 {
        self.values[row * self.size + column]
    }
}

struct SimilarityMatrix {
    columns: usize,
    values: Vec<f64>,
}

impl SimilarityMatrix {
    fn get(&self, row: usize, column: usize) -> f64 {
        self.values[row * self.columns + column]
    }
}

fn similarity_matrix(
    reference: &DistanceMatrix,
    mobile: &DistanceMatrix,
    window: usize,
) -> SimilarityMatrix {
    let rows = reference.size - window + 1;
    let columns = mobile.size - window + 1;
    let mut values = Vec::with_capacity(rows * columns);
    let terms = (window - 1) * (window - 2) / 2;
    for reference_start in 0..rows {
        for mobile_start in 0..columns {
            let mut difference = 0.0;
            for left in 0..window - 2 {
                for right in left + 2..window {
                    difference += (reference.get(reference_start + left, reference_start + right)
                        - mobile.get(mobile_start + left, mobile_start + right))
                    .abs();
                }
            }
            values.push(-difference / usize_to_f64(terms));
        }
    }
    SimilarityMatrix { columns, values }
}

#[derive(Clone, Copy)]
struct Fragment {
    reference: usize,
    mobile: usize,
}

struct Candidate {
    fragments: Vec<Fragment>,
    similarity: f64,
}

fn find_paths(
    similarity: &SimilarityMatrix,
    reference: &DistanceMatrix,
    mobile: &DistanceMatrix,
    options: CeOptions,
) -> Vec<Candidate> {
    let reference_starts = reference.size - options.window_size + 1;
    let mobile_starts = mobile.size - options.window_size + 1;
    let mut best = Vec::new();
    for reference_start in 0..reference_starts {
        if cannot_beat(reference_start, reference.size, &best, options.window_size) {
            break;
        }
        for mobile_start in 0..mobile_starts {
            if cannot_beat(mobile_start, mobile.size, &best, options.window_size) {
                break;
            }
            let initial = similarity.get(reference_start, mobile_start);
            if initial <= options.fragment_similarity_threshold {
                continue;
            }
            let mut candidate = Candidate {
                fragments: vec![Fragment {
                    reference: reference_start,
                    mobile: mobile_start,
                }],
                similarity: initial,
            };
            extend_path(&mut candidate, similarity, reference, mobile, options);
            insert_candidate(&mut best, candidate, options.max_paths);
        }
    }
    best
}

fn cannot_beat(start: usize, length: usize, best: &[Candidate], window: usize) -> bool {
    best.last().is_some_and(|worst| {
        let gaps = if worst.fragments.is_empty() {
            0
        } else {
            worst.fragments.len() - 1
        };
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
    candidate: &mut Candidate,
    similarity: &SimilarityMatrix,
    reference: &DistanceMatrix,
    mobile: &DistanceMatrix,
    options: CeOptions,
) {
    loop {
        let Some(last) = candidate.fragments.last().copied() else {
            return;
        };
        let mut best_extension: Option<(Fragment, f64)> = None;
        for gap_index in 0..=options.max_gap * 2 {
            let mut next = Fragment {
                reference: last.reference + options.window_size,
                mobile: last.mobile + options.window_size,
            };
            if (gap_index + 1).is_multiple_of(2) {
                next.reference += gap_index.div_ceil(2);
            } else {
                next.mobile += gap_index.div_ceil(2);
            }
            if next.reference + options.window_size > reference.size
                || next.mobile + options.window_size > mobile.size
                || similarity.get(next.reference, next.mobile)
                    <= options.fragment_similarity_threshold
            {
                continue;
            }
            let cross = candidate
                .fragments
                .iter()
                .map(|fragment| {
                    cross_similarity(*fragment, next, reference, mobile, options.window_size)
                })
                .sum::<f64>()
                / usize_to_f64(candidate.fragments.len());
            if cross > options.path_similarity_threshold
                && best_extension.is_none_or(|(_, best_cross)| cross > best_cross)
            {
                best_extension = Some((next, cross));
            }
        }
        let Some((next, cross)) = best_extension else {
            break;
        };
        let count = usize_to_f64(candidate.fragments.len());
        let current_terms = count + count * (count - 1.0) / 2.0;
        let new_terms = count + 1.0 + count * (count + 1.0) / 2.0;
        let next_similarity = (current_terms * candidate.similarity
            + count * cross
            + similarity.get(next.reference, next.mobile))
            / new_terms;
        if next_similarity <= options.path_similarity_threshold {
            break;
        }
        candidate.fragments.push(next);
        candidate.similarity = next_similarity;
    }
}

fn cross_similarity(
    first: Fragment,
    second: Fragment,
    reference: &DistanceMatrix,
    mobile: &DistanceMatrix,
    window: usize,
) -> f64 {
    let mut difference = (reference.get(first.reference, second.reference)
        - mobile.get(first.mobile, second.mobile))
    .abs();
    difference += (reference.get(first.reference + window - 1, second.reference + window - 1)
        - mobile.get(first.mobile + window - 1, second.mobile + window - 1))
    .abs();
    for offset in 1..window - 1 {
        difference += (reference.get(
            first.reference + offset,
            second.reference + window - 1 - offset,
        ) - mobile.get(first.mobile + offset, second.mobile + window - 1 - offset))
        .abs();
    }
    -difference / usize_to_f64(window)
}

fn insert_candidate(best: &mut Vec<Candidate>, candidate: Candidate, limit: usize) {
    best.push(candidate);
    best.sort_by(|left, right| {
        right
            .fragments
            .len()
            .cmp(&left.fragments.len())
            .then_with(|| right.similarity.total_cmp(&left.similarity))
    });
    best.truncate(limit);
}

fn finish(
    candidate: &Candidate,
    reference: &[[f32; 3]],
    mobile: &[[f32; 3]],
    options: CeOptions,
) -> Result<CeAlignment, CeError> {
    let window = options.window_size;
    let mut reference_indices = Vec::with_capacity(candidate.fragments.len() * window);
    let mut mobile_indices = Vec::with_capacity(candidate.fragments.len() * window);
    for fragment in &candidate.fragments {
        reference_indices.extend(fragment.reference..fragment.reference + window);
        mobile_indices.extend(fragment.mobile..fragment.mobile + window);
    }
    let reference_points: Vec<_> = reference_indices
        .iter()
        .map(|index| reference[*index])
        .collect();
    let mobile_points: Vec<_> = mobile_indices.iter().map(|index| mobile[*index]).collect();
    let rmsd = pdbiox_geom::superpose(&mobile_points, &reference_points)
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

fn significance(candidate: &Candidate, profile: CeSignificanceProfile) -> f64 {
    match profile {
        CeSignificanceProfile::OriginalWindowEight => original_window_eight_z_score(candidate),
    }
}

fn original_window_eight_z_score(candidate: &Candidate) -> f64 {
    const SIMILARITY_AVERAGES: [f64; 20] = [
        2.54, 2.51, 2.72, 3.01, 3.31, 3.61, 3.90, 4.19, 4.47, 4.74, 4.99, 5.22, 5.46, 5.70, 5.94,
        6.13, 6.36, 6.52, 6.68, 6.91,
    ];
    const SIMILARITY_SDS: [f64; 20] = [
        1.33, 0.88, 0.73, 0.71, 0.74, 0.80, 0.86, 0.92, 0.98, 1.04, 1.08, 1.10, 1.15, 1.19, 1.23,
        1.25, 1.32, 1.34, 1.36, 1.45,
    ];
    let length = candidate.fragments.len();
    let (average, deviation) = if length <= 20 {
        (SIMILARITY_AVERAGES[length - 1], SIMILARITY_SDS[length - 1])
    } else {
        (
            0.209_874 * usize_to_f64(length) + 2.944_714,
            0.039_487 * usize_to_f64(length) + 0.675_735,
        )
    };
    (average + candidate.similarity) / deviation
}

#[cfg(test)]
#[path = "ce_tests.rs"]
mod tests;
