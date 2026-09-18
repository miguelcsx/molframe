//! Public entry points over the shared affine dynamic program.

use super::dynamic::{Mode, run};
use super::linear::global_linear;
use super::types::{
    AlignError, Alignment, AlignmentMode, RegionAlignError, RegionError, RegionOptions,
};
use crate::matrix::SubstitutionMatrix;
use crate::scoring::Scoring;
use std::ops::Range;

/// Aligns two sequences end to end (Needleman–Wunsch with affine gaps).
///
/// Runs in `O(n·m)` time and space.
///
/// # Errors
///
/// Returns [`AlignError::NumericOverflow`] when dimensions or the exact score
/// exceed their supported representation.
pub fn global(left: &[u8], right: &[u8], scoring: Scoring) -> Result<Alignment, AlignError> {
    global_linear(left, right, &scoring, scoring.gap_open, scoring.gap_extend)
}

/// Aligns two sequences globally within a diagonal band of the given width.
///
/// Only cells whose row and column indices differ by at most `band` are
/// considered, dropping the cost to `O(n · band)`. The band must be at least the
/// length difference of the two sequences for a global path to exist; a narrower
/// band cannot reach the far corner and returns an explicit error.
///
/// A band wide enough to cover the whole matrix gives exactly [`global`].
///
/// # Errors
///
/// Returns [`AlignError::NumericOverflow`] when dimensions or the exact score
/// exceed their supported representation, and [`AlignError::NoAlignmentPath`]
/// when the band excludes every global path.
pub fn global_banded(
    left: &[u8],
    right: &[u8],
    scoring: Scoring,
    band: usize,
) -> Result<Alignment, AlignError> {
    run(
        left,
        right,
        &scoring,
        scoring.gap_open,
        scoring.gap_extend,
        Mode::Global,
        Some(band),
    )
}

/// Finds the highest-scoring local alignment (Smith–Waterman with affine gaps).
///
/// Runs in `O(n·m)` time and space.
///
/// # Errors
///
/// Returns [`AlignError::NumericOverflow`] when dimensions or the exact score
/// exceed their supported representation.
pub fn local(left: &[u8], right: &[u8], scoring: Scoring) -> Result<Alignment, AlignError> {
    run(
        left,
        right,
        &scoring,
        scoring.gap_open,
        scoring.gap_extend,
        Mode::Local,
        None,
    )
}

/// Aligns two sequences with the end gaps made free.
///
/// A short sequence aligned against a long one pays nothing for the unaligned
/// ends, so it settles onto its best internal match. Runs in `O(n·m)` time and
/// space.
///
/// # Errors
///
/// Returns [`AlignError::NumericOverflow`] when dimensions or the exact score
/// exceed their supported representation.
pub fn semi_global(left: &[u8], right: &[u8], scoring: Scoring) -> Result<Alignment, AlignError> {
    run(
        left,
        right,
        &scoring,
        scoring.gap_open,
        scoring.gap_extend,
        Mode::SemiGlobal,
        None,
    )
}

/// Aligns two explicit subranges and reports indices in the original sequences.
///
/// # Errors
///
/// Returns [`RegionAlignError`] when a range is invalid or the exact score
/// exceeds the supported numeric domain.
pub fn align_region(
    left: &[u8],
    right: &[u8],
    scoring: Scoring,
    options: &RegionOptions,
) -> Result<Alignment, RegionAlignError> {
    let left_region = region(left, &options.left, true)?;
    let right_region = region(right, &options.right, false)?;
    let mode = match options.mode {
        AlignmentMode::Global => Mode::Global,
        AlignmentMode::Local => Mode::Local,
        AlignmentMode::SemiGlobal => Mode::SemiGlobal,
    };
    let mut alignment = run(
        left_region,
        right_region,
        &scoring,
        scoring.gap_open,
        scoring.gap_extend,
        mode,
        options.band,
    )?;
    for column in &mut alignment.columns {
        column.left = column.left.map(|index| index + options.left.start);
        column.right = column.right.map(|index| index + options.right.start);
    }
    Ok(alignment)
}

fn region<'a>(
    sequence: &'a [u8],
    range: &Range<usize>,
    left: bool,
) -> Result<&'a [u8], RegionError> {
    sequence.get(range.clone()).ok_or(RegionError {
        left,
        length: sequence.len(),
    })
}

/// Aligns two protein sequences end to end scored by a substitution matrix.
///
/// The gap costs are given directly — `gap_open` is the score for the first
/// residue of a gap and `gap_extend` for each residue after it — because a
/// substitution matrix carries only substitution scores, not gap policy. Runs
/// in `O(n·m)` time and space.
///
/// # Examples
///
/// ```
/// use molframe_seq::align::global_matrix;
/// use molframe_seq::matrix::blosum62;
/// let alignment = global_matrix(b"WKFL", b"WKFL", &blosum62(), -11, -1)?;
/// // Two identical proteins align with no gaps and a positive score.
/// assert_eq!(alignment.columns.len(), 4);
/// assert!(alignment.score > 0);
/// # Ok::<(), molframe_seq::AlignError>(())
/// ```
///
/// # Errors
///
/// Returns [`AlignError::NumericOverflow`] when dimensions or the exact score
/// exceed their supported representation.
pub fn global_matrix(
    left: &[u8],
    right: &[u8],
    matrix: &SubstitutionMatrix,
    gap_open: i32,
    gap_extend: i32,
) -> Result<Alignment, AlignError> {
    global_linear(left, right, matrix, gap_open, gap_extend)
}

/// Finds the highest-scoring local alignment scored by a substitution matrix.
///
/// The counterpart of [`local`] for protein scoring. Runs in `O(n·m)` time and
/// space.
///
/// # Examples
///
/// ```
/// use molframe_seq::align::local_matrix;
/// use molframe_seq::matrix::blosum62;
/// let alignment = local_matrix(b"AAWKFAA", b"WKF", &blosum62(), -11, -1)?;
/// assert!(alignment.score > 0);
/// # Ok::<(), molframe_seq::AlignError>(())
/// ```
///
/// # Errors
///
/// Returns [`AlignError::NumericOverflow`] when dimensions or the exact score
/// exceed their supported representation.
pub fn local_matrix(
    left: &[u8],
    right: &[u8],
    matrix: &SubstitutionMatrix,
    gap_open: i32,
    gap_extend: i32,
) -> Result<Alignment, AlignError> {
    run(left, right, matrix, gap_open, gap_extend, Mode::Local, None)
}

/// Aligns two sequences with free end gaps, scored by a substitution matrix.
///
/// The counterpart of [`semi_global`] for protein scoring. Runs in `O(n·m)`
/// time and space.
///
/// # Examples
///
/// ```
/// use molframe_seq::align::semi_global_matrix;
/// use molframe_seq::matrix::blosum62;
/// let alignment = semi_global_matrix(b"GGWKFGG", b"WKF", &blosum62(), -11, -1)?;
/// assert!(alignment.score > 0);
/// # Ok::<(), molframe_seq::AlignError>(())
/// ```
///
/// # Errors
///
/// Returns [`AlignError::NumericOverflow`] when dimensions or the exact score
/// exceed their supported representation.
pub fn semi_global_matrix(
    left: &[u8],
    right: &[u8],
    matrix: &SubstitutionMatrix,
    gap_open: i32,
    gap_extend: i32,
) -> Result<Alignment, AlignError> {
    run(
        left,
        right,
        matrix,
        gap_open,
        gap_extend,
        Mode::SemiGlobal,
        None,
    )
}

/// Scores a global affine alignment without constructing its traceback.
///
/// Uses `O(min(n, m))` working memory and `O(n * m)` time.
///
/// # Errors
///
/// Returns a numeric error if dimensions or the exact score overflow.
pub fn global_score(left: &[u8], right: &[u8], scoring: Scoring) -> Result<i32, AlignError> {
    super::linear::global_score_linear(left, right, &scoring, scoring.gap_open, scoring.gap_extend)
}
