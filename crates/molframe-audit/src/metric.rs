//! How far apart two results of the same analysis are.
//!
//! An audit varies decisions and asks how much the answer moved. "How much"
//! depends on what the answer is: a set of contacts, one number, a profile, a
//! ranking, a network. Each kind has its own natural distance, and the audit
//! needs only that distance, so the kinds share one trait and the audit never
//! learns which one it is looking at.
//!
//! Every metric returns a non-negative, symmetric number that is zero for
//! identical results. Where two results cannot be compared at all (a vector
//! whose length changed, a number that became `NaN`) the distance is infinite,
//! which is the honest value: an analysis whose output changes shape under a
//! decision is not stable under it, and a finite stand-in would hide that.

use crate::numeric::usize_to_f64;
use std::collections::{BTreeMap, BTreeSet};

/// A distance between two results of one analysis.
pub trait OutcomeMetric<R> {
    /// A short stable name, recorded with every impact this metric measures.
    fn name(&self) -> &'static str;

    /// The distance between two results; symmetric, non-negative, zero when equal.
    fn distance(&self, first: &R, second: &R) -> f64;
}

/// One minus the Jaccard index of two sets of identities.
///
/// For contacts, hydrogen bonds, residues, pockets and any other result that is
/// "which items are present". Two empty sets are identical.
pub struct SetOverlap<P> {
    items: P,
}

impl<P> std::fmt::Debug for SetOverlap<P> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.debug_struct("SetOverlap").finish_non_exhaustive()
    }
}

impl<P> SetOverlap<P> {
    /// Compares the items `items` projects out of each result.
    pub const fn new(items: P) -> Self {
        Self { items }
    }
}

impl<R, I: Ord, P: Fn(&R) -> BTreeSet<I>> OutcomeMetric<R> for SetOverlap<P> {
    fn name(&self) -> &'static str {
        "jaccard-distance"
    }

    fn distance(&self, first: &R, second: &R) -> f64 {
        jaccard_distance(&(self.items)(first), &(self.items)(second))
    }
}

/// One minus the Jaccard index of two sets.
pub(crate) fn jaccard_distance<I: Ord>(first: &BTreeSet<I>, second: &BTreeSet<I>) -> f64 {
    let union = first.union(second).count();
    if union == 0 {
        0.0
    } else {
        1.0 - usize_to_f64(first.intersection(second).count()) / usize_to_f64(union)
    }
}

/// Whether the conclusion changed: zero when two results fall in the same category, one
/// when they do not.
///
/// A score that moves from 0.64 to 0.61 has changed little; the same score crossing the
/// line between "acceptable" and "incorrect" has changed its conclusion. Projecting each
/// result to the category a threshold puts it in measures that, and the mean distance over
/// pairs of runs is the probability that two defensible universes disagree.
pub struct CategoricalFlip<P> {
    category: P,
}

impl<P> std::fmt::Debug for CategoricalFlip<P> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("CategoricalFlip")
            .finish_non_exhaustive()
    }
}

impl<P> CategoricalFlip<P> {
    /// Compares the category `category` projects out of each result.
    pub const fn new(category: P) -> Self {
        Self { category }
    }
}

impl<R, C: PartialEq, P: Fn(&R) -> C> OutcomeMetric<R> for CategoricalFlip<P> {
    fn name(&self) -> &'static str {
        "conclusion-flip"
    }

    fn distance(&self, first: &R, second: &R) -> f64 {
        if (self.category)(first) == (self.category)(second) {
            0.0
        } else {
            1.0
        }
    }
}

/// How a difference between two numbers is measured.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ScalarMode {
    /// `|a - b|`, in the units of the quantity.
    Absolute,
    /// `|a - b| / max(|a|, |b|)`, unitless; zero when both are zero.
    Relative,
}

/// The difference between two numbers.
pub struct ScalarError<P> {
    value: P,
    mode: ScalarMode,
}

impl<P> std::fmt::Debug for ScalarError<P> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("ScalarError")
            .finish_non_exhaustive()
    }
}

impl<P> ScalarError<P> {
    /// Compares the number `value` projects out of each result.
    pub const fn new(value: P, mode: ScalarMode) -> Self {
        Self { value, mode }
    }
}

impl<R, P: Fn(&R) -> f64> OutcomeMetric<R> for ScalarError<P> {
    fn name(&self) -> &'static str {
        match self.mode {
            ScalarMode::Absolute => "absolute-error",
            ScalarMode::Relative => "relative-error",
        }
    }

    fn distance(&self, first: &R, second: &R) -> f64 {
        let (a, b) = ((self.value)(first), (self.value)(second));
        if !a.is_finite() || !b.is_finite() {
            // Two identical non-finite values agree; anything else cannot be compared.
            return if a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan()) {
                0.0
            } else {
                f64::INFINITY
            };
        }
        let difference = (a - b).abs();
        match self.mode {
            ScalarMode::Absolute => difference,
            ScalarMode::Relative => {
                let scale = a.abs().max(b.abs());
                if scale == 0.0 {
                    0.0
                } else {
                    difference / scale
                }
            }
        }
    }
}

/// How two vectors of the same length are compared.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum VectorMode {
    /// Root-mean-square difference, in the units of the quantity.
    Rms,
    /// `sqrt(2 (1 - r))` for Pearson's `r`: the Euclidean distance between the
    /// standardised vectors, so a metric where `1 - r` is not.
    ///
    /// Two constant vectors agree; a constant against a varying one is as far
    /// as uncorrelated vectors are.
    Correlation,
}

/// The difference between two vectors, such as a per-residue profile.
pub struct VectorDifference<P> {
    values: P,
    mode: VectorMode,
}

impl<P> std::fmt::Debug for VectorDifference<P> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("VectorDifference")
            .finish_non_exhaustive()
    }
}

impl<P> VectorDifference<P> {
    /// Compares the vector `values` projects out of each result.
    pub const fn new(values: P, mode: VectorMode) -> Self {
        Self { values, mode }
    }
}

impl<R, P: Fn(&R) -> Vec<f64>> OutcomeMetric<R> for VectorDifference<P> {
    fn name(&self) -> &'static str {
        match self.mode {
            VectorMode::Rms => "rms-difference",
            VectorMode::Correlation => "correlation-distance",
        }
    }

    fn distance(&self, first: &R, second: &R) -> f64 {
        let (a, b) = ((self.values)(first), (self.values)(second));
        if a.len() != b.len() || a.iter().chain(&b).any(|value| !value.is_finite()) {
            return f64::INFINITY;
        }
        if a.is_empty() {
            return 0.0;
        }
        match self.mode {
            VectorMode::Rms => {
                let squares: f64 = a.iter().zip(&b).map(|(x, y)| (x - y).powi(2)).sum();
                (squares / usize_to_f64(a.len())).sqrt()
            }
            VectorMode::Correlation => correlation_distance(&a, &b),
        }
    }
}

fn correlation_distance(a: &[f64], b: &[f64]) -> f64 {
    let n = usize_to_f64(a.len());
    let (mean_a, mean_b) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
    let (mut covariance, mut variance_a, mut variance_b) = (0.0, 0.0, 0.0);
    for (x, y) in a.iter().zip(b) {
        covariance += (x - mean_a) * (y - mean_b);
        variance_a += (x - mean_a).powi(2);
        variance_b += (y - mean_b).powi(2);
    }
    match (variance_a == 0.0, variance_b == 0.0) {
        (true, true) => 0.0,
        (true, false) | (false, true) => 2.0_f64.sqrt(),
        (false, false) => {
            let r = (covariance / (variance_a * variance_b).sqrt()).clamp(-1.0, 1.0);
            (2.0 * (1.0 - r)).max(0.0).sqrt()
        }
    }
}

/// The disagreement between two rankings of items, best first.
///
/// Kendall's tau distance for partial rankings (Fagin et al., 2006) with the
/// neutral penalty `p = 1/2`: an item a ranking does not list is tied, below
/// every item it does list, with the other unlisted ones. A pair ordered
/// oppositely costs one, a pair ordered by one ranking and tied in the other
/// costs one half, and the total is divided by the number of pairs, so the
/// distance lies in `[0, 1]`.
pub struct RankingDistance<P> {
    ranking: P,
}

impl<P> std::fmt::Debug for RankingDistance<P> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("RankingDistance")
            .finish_non_exhaustive()
    }
}

impl<P> RankingDistance<P> {
    /// Compares the best-first list `ranking` projects out of each result.
    pub const fn new(ranking: P) -> Self {
        Self { ranking }
    }
}

impl<R, I: Ord + Clone, P: Fn(&R) -> Vec<I>> OutcomeMetric<R> for RankingDistance<P> {
    fn name(&self) -> &'static str {
        "kendall-tau-distance"
    }

    fn distance(&self, first: &R, second: &R) -> f64 {
        let (a, b) = ((self.ranking)(first), (self.ranking)(second));
        kendall_distance(&positions(&a), &positions(&b))
    }
}

fn positions<I: Ord + Clone>(ranking: &[I]) -> BTreeMap<I, usize> {
    let mut positions = BTreeMap::new();
    for (rank, item) in ranking.iter().enumerate() {
        positions.entry(item.clone()).or_insert(rank);
    }
    positions
}

fn kendall_distance<I: Ord + Clone>(a: &BTreeMap<I, usize>, b: &BTreeMap<I, usize>) -> f64 {
    let items: Vec<&I> = a
        .keys()
        .chain(b.keys())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let mut penalty = 0.0;
    let mut pairs = 0.0;
    for (index, &x) in items.iter().enumerate() {
        for &y in &items[index + 1..] {
            pairs += 1.0;
            penalty += pair_penalty(order(a, x, y), order(b, x, y));
        }
    }
    if pairs == 0.0 { 0.0 } else { penalty / pairs }
}

/// How one ranking orders `x` against `y`: `Some(true)` when `x` comes first,
/// `Some(false)` when `y` does, `None` when it ties them (both unlisted).
fn order<I: Ord>(ranking: &BTreeMap<I, usize>, x: &I, y: &I) -> Option<bool> {
    match (ranking.get(x), ranking.get(y)) {
        (Some(rx), Some(ry)) => Some(rx < ry),
        (Some(_), None) => Some(true),
        (None, Some(_)) => Some(false),
        (None, None) => None,
    }
}

fn pair_penalty(first: Option<bool>, second: Option<bool>) -> f64 {
    match (first, second) {
        (Some(a), Some(b)) => {
            if a == b {
                0.0
            } else {
                1.0
            }
        }
        (None, None) => 0.0,
        (Some(_), None) | (None, Some(_)) => 0.5,
    }
}

/// A network as the identities of its nodes and of its edges.
///
/// Edges are unordered pairs: `(a, b)` and `(b, a)` are the same edge.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Graph<N: Ord> {
    nodes: BTreeSet<N>,
    edges: BTreeSet<(N, N)>,
}

impl<N: Ord + Clone> Graph<N> {
    /// An empty network.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            nodes: BTreeSet::new(),
            edges: BTreeSet::new(),
        }
    }

    /// Adds a node.
    pub fn add_node(&mut self, node: N) {
        self.nodes.insert(node);
    }

    /// Adds an edge, and its two endpoints as nodes.
    pub fn add_edge(&mut self, first: N, second: N) {
        self.nodes.insert(first.clone());
        self.nodes.insert(second.clone());
        let edge = if first <= second {
            (first, second)
        } else {
            (second, first)
        };
        self.edges.insert(edge);
    }

    /// The nodes.
    #[must_use]
    pub const fn nodes(&self) -> &BTreeSet<N> {
        &self.nodes
    }

    /// The edges, each as `(smaller, larger)`.
    #[must_use]
    pub const fn edges(&self) -> &BTreeSet<(N, N)> {
        &self.edges
    }

    /// One minus the Jaccard index of the two node sets.
    #[must_use]
    pub fn node_difference(&self, other: &Self) -> f64 {
        jaccard_distance(&self.nodes, &other.nodes)
    }

    /// One minus the Jaccard index of the two edge sets.
    #[must_use]
    pub fn edge_difference(&self, other: &Self) -> f64 {
        jaccard_distance(&self.edges, &other.edges)
    }
}

/// Which part of a network is compared.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum GraphPart {
    /// Which nodes exist.
    Nodes,
    /// Which edges exist.
    Edges,
    /// Nodes and edges together, as one set of identities.
    Both,
}

/// The difference between two networks.
pub struct GraphDifference<P> {
    graph: P,
    part: GraphPart,
}

impl<P> std::fmt::Debug for GraphDifference<P> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("GraphDifference")
            .finish_non_exhaustive()
    }
}

impl<P> GraphDifference<P> {
    /// Compares the `part` of the network `graph` projects out of each result.
    pub const fn new(graph: P, part: GraphPart) -> Self {
        Self { graph, part }
    }
}

impl<R, N: Ord + Clone, P: Fn(&R) -> Graph<N>> OutcomeMetric<R> for GraphDifference<P> {
    fn name(&self) -> &'static str {
        match self.part {
            GraphPart::Nodes => "node-difference",
            GraphPart::Edges => "edge-difference",
            GraphPart::Both => "graph-difference",
        }
    }

    fn distance(&self, first: &R, second: &R) -> f64 {
        let (a, b) = ((self.graph)(first), (self.graph)(second));
        match self.part {
            GraphPart::Nodes => a.node_difference(&b),
            GraphPart::Edges => a.edge_difference(&b),
            GraphPart::Both => {
                let identities = |graph: &Graph<N>| -> BTreeSet<(bool, N, Option<N>)> {
                    graph
                        .nodes
                        .iter()
                        .map(|node| (false, node.clone(), None))
                        .chain(
                            graph
                                .edges
                                .iter()
                                .map(|(x, y)| (true, x.clone(), Some(y.clone()))),
                        )
                        .collect()
                };
                jaccard_distance(&identities(&a), &identities(&b))
            }
        }
    }
}

#[cfg(test)]
#[path = "metric_tests.rs"]
mod tests;
