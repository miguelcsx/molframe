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

/// Builds a neighbour-joining tree from labels and their distance matrix.
///
/// Returns `None` when there are no taxa or the matrix is not square with the
/// labels. One taxon yields a lone leaf; two are joined directly.
///
/// Runs in `O(n³)` time.
#[must_use]
pub fn neighbor_joining(labels: &[&str], distances: &[Vec<f64>]) -> Option<Tree> {
    let n = labels.len();
    if n == 0 || distances.len() != n || distances.iter().any(|row| row.len() != n) {
        return None;
    }
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
