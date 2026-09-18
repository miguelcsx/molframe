//! UPGMA clustering of a distance matrix into a tree.
//!
//! UPGMA repeatedly joins the two closest clusters, placing the new node at half
//! their distance and averaging distances to the merged cluster by size. It
//! assumes a roughly constant rate of change, which makes it a simple, fully
//! deterministic way to turn a table of pairwise distances into a rooted tree.
//!
//! Ties are broken toward the lowest-indexed pair, so the same matrix always
//! yields the same tree. Cost is `O(n³)` in the number of taxa.

use crate::tree::Tree;

/// One live cluster: its subtree, the number of taxa in it, and its height.
struct Cluster {
    tree: Tree,
    size: f64,
    height: f64,
}

/// Builds a UPGMA tree from taxon labels and their symmetric distance matrix.
///
/// Returns `None` when there are no taxa or the matrix is not square with the
/// labels. A single taxon yields a lone leaf.
///
/// Runs in `O(n³)` time.
#[must_use]
pub fn upgma(labels: &[&str], distances: &[Vec<f64>]) -> Option<Tree> {
    let n = labels.len();
    if n == 0 || distances.len() != n || distances.iter().any(|row| row.len() != n) {
        return None;
    }

    let mut clusters: Vec<Option<Cluster>> = labels
        .iter()
        .map(|label| {
            Some(Cluster {
                tree: Tree::Leaf {
                    name: (*label).to_string(),
                },
                size: 1.0,
                height: 0.0,
            })
        })
        .collect();
    let mut distance = distances.to_vec();
    let mut alive = vec![true; n];
    let mut remaining = n;

    while remaining > 1 {
        let Some((_, i, j)) = closest_pair(&distance, &alive) else {
            break;
        };
        let height = distance[i][j] / 2.0;
        let left = clusters[i].take()?;
        let right = clusters[j].take()?;
        let merged = Cluster {
            tree: Tree::Clade {
                left: Box::new(left.tree),
                left_length: height - left.height,
                right: Box::new(right.tree),
                right_length: height - right.height,
            },
            size: left.size + right.size,
            height,
        };
        for k in 0..n {
            if !alive[k] || k == i || k == j {
                continue;
            }
            let averaged = (left.size * distance[i][k] + right.size * distance[j][k]) / merged.size;
            distance[i][k] = averaged;
            distance[k][i] = averaged;
        }
        clusters[i] = Some(merged);
        alive[j] = false;
        remaining -= 1;
    }

    let surviving = alive.iter().position(|&live| live)?;
    clusters[surviving].take().map(|cluster| cluster.tree)
}

/// The closest live pair `(distance, i, j)` with `i < j`, lowest indices first.
fn closest_pair(distance: &[Vec<f64>], alive: &[bool]) -> Option<(f64, usize, usize)> {
    let mut best: Option<(f64, usize, usize)> = None;
    for i in 0..alive.len() {
        if !alive[i] {
            continue;
        }
        for j in (i + 1)..alive.len() {
            if !alive[j] {
                continue;
            }
            let candidate = distance[i][j];
            match best {
                Some((current, _, _)) if candidate >= current => {}
                _ => best = Some((candidate, i, j)),
            }
        }
    }
    best
}

#[cfg(test)]
#[path = "algorithm_tests.rs"]
mod tests;
