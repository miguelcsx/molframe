//! Combinatorial Extension (CE) structural alignment over guide coordinates.
//!
//! Fragment geometry is retained as narrow distance bands rather than complete
//! pairwise matrices. The `R x M` fragment-similarity table is cached only when
//! it fits beside the complete search workspace under the caller's ceiling.

#[path = "ce/finish.rs"]
mod finish;
#[path = "ce/search.rs"]
mod search;
#[path = "ce/workspace.rs"]
mod workspace;

use pdbiox_geom::SuperposeError;

use crate::numeric::usize_to_f64;

/// Default byte ceiling for one CE search workspace: 100 MB.
///
/// A default, not a maximum. A caller aligning very long chains may raise it.
const DEFAULT_CE_MEMORY_LIMIT_BYTES: usize = 100_000_000;

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
    /// Operation-owned byte ceiling for the search workspace.
    pub memory_limit_bytes: usize,
}

impl CeOptions {
    /// Original CE search and significance profile under the default ceiling.
    #[must_use]
    pub const fn original() -> Self {
        Self {
            window_size: 8,
            max_gap: 30,
            max_paths: 20,
            fragment_similarity_threshold: -3.0,
            path_similarity_threshold: -4.0,
            significance: Some(CeSignificanceProfile::OriginalWindowEight),
            memory_limit_bytes: DEFAULT_CE_MEMORY_LIMIT_BYTES,
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

/// Invalid CE settings, inputs, allocation or final fit.
#[derive(Clone, Copy, Debug, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum CeError {
    /// Window, gap, path count or memory ceiling is outside the supported domain.
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
    /// A search or fit phase cannot remain within the caller's memory ceiling.
    #[error(
        "CE cannot provide an accounted {required}-byte workspace within the {limit}-byte limit"
    )]
    MemoryLimit {
        /// Required operation-owned bytes, or `usize::MAX` after overflow.
        required: usize,
        /// Caller-provided ceiling.
        limit: usize,
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
/// Returns an option, input, memory, no-path or superposition error.
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

/// Returns retained CE candidates ordered by coverage and fitted RMSD.
///
/// Runtime is `O(R M W² + R M P (G + P W))` in guide atoms, retained path
/// length `P`, gap bound `G` and window `W`. Search memory is
/// `O((R + M) W + P max_paths)` plus an optional, budgeted `O(R M)` cache.
///
/// # Errors
///
/// Returns an option, input, memory, no-path or superposition error.
pub fn ce_alignments(
    reference: &[[f32; 3]],
    mobile: &[[f32; 3]],
    options: CeOptions,
) -> Result<Vec<CeAlignment>, CeError> {
    validate(reference, mobile, options)?;
    let candidates = search::find_paths(reference, mobile, options)?;
    if candidates.is_empty() {
        return Err(CeError::NoAlignment);
    }
    finish::finish_candidates(candidates, reference, mobile, options)
}

fn validate(
    reference: &[[f32; 3]],
    mobile: &[[f32; 3]],
    options: CeOptions,
) -> Result<(), CeError> {
    if options.window_size < 3
        || options.max_paths == 0
        || options.memory_limit_bytes == 0
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

fn significance(candidate: &search::Candidate, profile: CeSignificanceProfile) -> f64 {
    match profile {
        CeSignificanceProfile::OriginalWindowEight => original_window_eight_z_score(candidate),
    }
}

fn original_window_eight_z_score(candidate: &search::Candidate) -> f64 {
    const AVERAGES: [f64; 20] = [
        2.54, 2.51, 2.72, 3.01, 3.31, 3.61, 3.90, 4.19, 4.47, 4.74, 4.99, 5.22, 5.46, 5.70, 5.94,
        6.13, 6.36, 6.52, 6.68, 6.91,
    ];
    const DEVIATIONS: [f64; 20] = [
        1.33, 0.88, 0.73, 0.71, 0.74, 0.80, 0.86, 0.92, 0.98, 1.04, 1.08, 1.10, 1.15, 1.19, 1.23,
        1.25, 1.32, 1.34, 1.36, 1.45,
    ];
    let length = candidate.fragments.len();
    let (average, deviation) = if length <= 20 {
        (AVERAGES[length - 1], DEVIATIONS[length - 1])
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

#[cfg(test)]
#[path = "ce_reference_tests.rs"]
mod reference_tests;
