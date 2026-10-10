//! Neighbour-joining of a distance matrix into a tree.
//!
//! Neighbour joining builds a tree without assuming a constant rate of change:
//! at each step it joins the pair that most reduces the total branch length,
//! found through the Q-matrix, and assigns branch lengths from the row sums. On
//! additive distances — those that come from a real tree — it recovers that tree
//! exactly, which UPGMA cannot promise.
//!
//! The result is an unrooted tree; it is returned rooted at the last cluster
//! joined, with that cluster's branch given length zero. Ties are broken toward
//! the lowest-indexed pair, so the same matrix always yields the same tree. Cost
//! is `O(n³)`.

use crate::tree::Tree;

/// Absolute-plus-relative tolerance for matrix symmetry.
const SYMMETRY_TOLERANCE: f64 = 1.0e-9;

/// Why a distance matrix cannot be joined.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NjError {
    /// No taxa were given.
    Empty,
    /// The matrix is not square with the labels.
    NotSquare,
    /// An entry is NaN or infinite.
    NotFinite {
        /// Row of the offending entry.
        row: usize,
        /// Column of the offending entry.
        column: usize,
    },
    /// An entry is negative.
    Negative {
        /// Row of the offending entry.
        row: usize,
        /// Column of the offending entry.
        column: usize,
    },
    /// Entries `(row, column)` and `(column, row)` differ.
    Asymmetric {
        /// Row of the offending entry.
        row: usize,
        /// Column of the offending entry.
        column: usize,
    },
    /// A diagonal entry is not zero.
    NonZeroDiagonal {
        /// Index of the offending diagonal entry.
        index: usize,
    },
}

impl std::fmt::Display for NjError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Empty => write!(formatter, "neighbour joining needs at least one taxon"),
            Self::NotSquare => write!(
                formatter,
                "the distance matrix is not square with the labels"
            ),
            Self::NotFinite { row, column } => {
                write!(formatter, "distance ({row}, {column}) is not finite")
            }
            Self::Negative { row, column } => {
                write!(formatter, "distance ({row}, {column}) is negative")
            }
            Self::Asymmetric { row, column } => {
                write!(
                    formatter,
                    "distances ({row}, {column}) and ({column}, {row}) differ"
                )
            }
            Self::NonZeroDiagonal { index } => {
                write!(formatter, "diagonal distance {index} is not zero")
            }
        }
    }
}

impl std::error::Error for NjError {}

/// What to do when the algorithm computes a negative branch length.
///
/// Neighbour joining can yield negative lengths from perfectly valid but
/// non-additive distances. They are an artefact of the estimator, not a
/// meaningful quantity, yet clamping them changes the tree's path lengths.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum NegativeBranches {
    /// Report the computed value unchanged (the default, which keeps path
    /// lengths consistent with the algorithm's arithmetic).
    #[default]
    Keep,
    /// Replace negative lengths with zero.
    ClampToZero,
}

/// Options for [`neighbor_joining_with`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct NjOptions {
    /// Handling of negative branch lengths.
    pub negative_branches: NegativeBranches,
}

/// Checks that the matrix is square, finite, non-negative, symmetric within
/// tolerance and has a zero diagonal.
///
/// # Errors
///
/// Returns the first violated condition.
pub fn validate_distance_matrix(labels: &[&str], distances: &[Vec<f64>]) -> Result<(), NjError> {
    let n = labels.len();
    if n == 0 {
        return Err(NjError::Empty);
    }
    if distances.len() != n || distances.iter().any(|row| row.len() != n) {
        return Err(NjError::NotSquare);
    }
    for (row, values) in distances.iter().enumerate() {
        for (column, value) in values.iter().enumerate() {
            if !value.is_finite() {
                return Err(NjError::NotFinite { row, column });
            }
            if *value < 0.0 {
                return Err(NjError::Negative { row, column });
            }
        }
        if values[row] != 0.0 {
            return Err(NjError::NonZeroDiagonal { index: row });
        }
    }
    for (row, values) in distances.iter().enumerate() {
        for (column, left) in values.iter().enumerate().skip(row + 1) {
            let (left, right) = (*left, distances[column][row]);
            if (left - right).abs() > SYMMETRY_TOLERANCE * (1.0 + left.abs().max(right.abs())) {
                return Err(NjError::Asymmetric { row, column });
            }
        }
    }
    Ok(())
}

/// Builds a neighbour-joining tree with default options.
///
/// Returns `None` when the matrix fails [`validate_distance_matrix`]; use
/// [`neighbor_joining_with`] to learn why. One taxon yields a lone leaf; two
/// are joined directly. Negative branch lengths are kept.
///
/// Runs in `O(n³)` time.
#[must_use]
pub fn neighbor_joining(labels: &[&str], distances: &[Vec<f64>]) -> Option<Tree> {
    neighbor_joining_with(labels, distances, NjOptions::default()).ok()
}

/// Builds a neighbour-joining tree from labels and their distance matrix.
///
/// # Errors
///
/// Returns an [`NjError`] when the matrix is empty, ragged, non-finite,
/// negative, asymmetric beyond tolerance or has a non-zero diagonal.
///
/// Runs in `O(n³)` time.
pub fn neighbor_joining_with(
    labels: &[&str],
    distances: &[Vec<f64>],
    options: NjOptions,
) -> Result<Tree, NjError> {
    validate_distance_matrix(labels, distances)?;
    let tree = join(labels, distances).ok_or(NjError::NotSquare)?;
    Ok(match options.negative_branches {
        NegativeBranches::Keep => tree,
        NegativeBranches::ClampToZero => clamp_negative(tree),
    })
}

fn clamp_negative(tree: Tree) -> Tree {
    match tree {
        Tree::Leaf { .. } => tree,
        Tree::Clade {
            left,
            left_length,
            right,
            right_length,
        } => Tree::Clade {
            left: Box::new(clamp_negative(*left)),
            left_length: left_length.max(0.0),
            right: Box::new(clamp_negative(*right)),
            right_length: right_length.max(0.0),
        },
    }
}

fn join(labels: &[&str], distances: &[Vec<f64>]) -> Option<Tree> {
    let n = labels.len();
    if n == 1 {
        return Some(Tree::Leaf {
            name: labels[0].to_string(),
        });
    }

    let mut trees: Vec<Option<Tree>> = labels
        .iter()
        .map(|label| {
            Some(Tree::Leaf {
                name: (*label).to_string(),
            })
        })
        .collect();
    let mut distance = distances.to_vec();
    let mut alive = vec![true; n];
    let mut remaining = u32::try_from(n).ok()?;

    while remaining > 2 {
        let sums = row_sums(&distance, &alive);
        let (i, j) = lowest_q(&distance, &alive, &sums, remaining)?;
        let separation = distance[i][j];
        let delta = 0.5 * separation + (sums[i] - sums[j]) / (2.0 * (f64::from(remaining) - 2.0));
        let left = trees[i].take()?;
        let right = trees[j].take()?;
        trees[i] = Some(Tree::Clade {
            left: Box::new(left),
            left_length: delta,
            right: Box::new(right),
            right_length: separation - delta,
        });
        for k in 0..n {
            if !alive[k] || k == i || k == j {
                continue;
            }
            let updated = 0.5 * (distance[i][k] + distance[j][k] - separation);
            distance[i][k] = updated;
            distance[k][i] = updated;
        }
        alive[j] = false;
        remaining -= 1;
    }

    let (i, j) = live_pair(&alive)?;
    let left = trees[i].take()?;
    let right = trees[j].take()?;
    Some(Tree::Clade {
        left: Box::new(left),
        left_length: 0.0,
        right: Box::new(right),
        right_length: distance[i][j],
    })
}

/// Sum of distances from each live cluster to every other.
fn row_sums(distance: &[Vec<f64>], alive: &[bool]) -> Vec<f64> {
    let mut sums = vec![0.0; alive.len()];
    for i in 0..alive.len() {
        if !alive[i] {
            continue;
        }
        for k in 0..alive.len() {
            if alive[k] && k != i {
                sums[i] += distance[i][k];
            }
        }
    }
    sums
}

/// The live pair minimising the Q-criterion, lowest indices breaking ties.
fn lowest_q(
    distance: &[Vec<f64>],
    alive: &[bool],
    sums: &[f64],
    remaining: u32,
) -> Option<(usize, usize)> {
    let scale = f64::from(remaining) - 2.0;
    let mut best: Option<(f64, usize, usize)> = None;
    for i in 0..alive.len() {
        if !alive[i] {
            continue;
        }
        for j in (i + 1)..alive.len() {
            if !alive[j] {
                continue;
            }
            let q = scale * distance[i][j] - sums[i] - sums[j];
            match best {
                Some((current, _, _)) if q >= current => {}
                _ => best = Some((q, i, j)),
            }
        }
    }
    best.map(|(_, i, j)| (i, j))
}

/// The two live indices when exactly two remain.
fn live_pair(alive: &[bool]) -> Option<(usize, usize)> {
    let mut live = (0..alive.len()).filter(|&index| alive[index]);
    let first = live.next()?;
    let second = live.next()?;
    Some((first, second))
}

#[cfg(test)]
#[path = "algorithm_tests.rs"]
mod tests;
