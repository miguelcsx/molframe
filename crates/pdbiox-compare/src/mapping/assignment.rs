//! Maximum-total-identity one-to-one chain assignment.

use super::chain::{ChainSequence, chain_sequences, identity};
use pdbiox_chem::ComponentProvider;
use pdbiox_core::contract::Namespace;
use pdbiox_core::structure::Structure;
use pdbiox_core::{Code, Diagnostic};
use pdbiox_seq::Scoring;

/// A matched pair of chains and how identical their sequences are.
#[derive(Clone, Debug, PartialEq)]
pub struct ChainMapping {
    /// The reference chain label.
    pub reference: String,
    /// The target chain label.
    pub target: String,
    /// Fraction of aligned positions that are identical, from 0 to 1.
    pub identity: f64,
}

/// A viable non-primary target for one reference chain.
#[derive(Clone, Debug, PartialEq)]
pub struct ChainAlternative {
    /// The reference chain label.
    pub reference: String,
    /// An alternative target chain label.
    pub target: String,
    /// Pairwise aligned sequence identity.
    pub identity: f64,
}

/// Globally optimal chain assignment plus every viable pairwise alternative.
#[derive(Clone, Debug, PartialEq)]
pub struct ChainAssignment {
    /// Maximum-total-identity one-to-one assignment.
    pub primary: Vec<ChainMapping>,
    /// Unselected pairs meeting the same explicit identity threshold.
    pub alternatives: Vec<ChainAlternative>,
}

/// Maps reference chains to targets with a globally optimal assignment.
///
/// # Errors
///
/// Returns a CCD, namespace, or identity-threshold diagnostic.
pub fn map_chains(
    reference: &Structure,
    target: &Structure,
    provider: &dyn ComponentProvider,
    namespace: Namespace,
    scoring: Scoring,
    min_identity: f64,
) -> Result<Vec<ChainMapping>, Diagnostic> {
    assign_chains(
        reference,
        target,
        provider,
        namespace,
        scoring,
        min_identity,
    )
    .map(|assignment| assignment.primary)
}

/// Returns the optimal assignment and all threshold-qualified pair alternatives.
///
/// # Errors
///
/// Returns a CCD, namespace, or identity-threshold diagnostic.
pub fn assign_chains(
    reference: &Structure,
    target: &Structure,
    provider: &dyn ComponentProvider,
    namespace: Namespace,
    scoring: Scoring,
    min_identity: f64,
) -> Result<ChainAssignment, Diagnostic> {
    validate_threshold(min_identity)?;
    let references = chain_sequences(reference, provider, namespace)?;
    let targets = chain_sequences(target, provider, namespace)?;
    let scores = score_matrix(&references, &targets, scoring)?;
    let selected = optimal_pairs(&scores, min_identity);
    Ok(ChainAssignment {
        primary: primary_mappings(&references, &targets, &scores, &selected),
        alternatives: alternatives(&references, &targets, &scores, &selected, min_identity),
    })
}

fn validate_threshold(value: f64) -> Result<(), Diagnostic> {
    if value.is_finite() && (0.0..=1.0).contains(&value) {
        Ok(())
    } else {
        Err(Diagnostic::new(Code::E4003).with_context("min_identity", value.to_string()))
    }
}

fn score_matrix(
    references: &[ChainSequence],
    targets: &[ChainSequence],
    scoring: Scoring,
) -> Result<Vec<Vec<f64>>, Diagnostic> {
    references
        .iter()
        .map(|reference| {
            targets
                .iter()
                .map(|target| identity(&reference.sequence, &target.sequence, scoring))
                .collect()
        })
        .collect()
}

fn primary_mappings(
    references: &[ChainSequence],
    targets: &[ChainSequence],
    scores: &[Vec<f64>],
    selected: &[(usize, usize)],
) -> Vec<ChainMapping> {
    let mut mappings: Vec<_> = selected
        .iter()
        .map(|&(reference, target)| ChainMapping {
            reference: references[reference].label.clone(),
            target: targets[target].label.clone(),
            identity: scores[reference][target],
        })
        .collect();
    mappings.sort_by(|left, right| left.reference.cmp(&right.reference));
    mappings
}

fn alternatives(
    references: &[ChainSequence],
    targets: &[ChainSequence],
    scores: &[Vec<f64>],
    selected: &[(usize, usize)],
    threshold: f64,
) -> Vec<ChainAlternative> {
    let mut alternatives = Vec::new();
    for (reference, row) in scores.iter().enumerate() {
        for (target, &pair_identity) in row.iter().enumerate() {
            if pair_identity >= threshold && !selected.contains(&(reference, target)) {
                alternatives.push(ChainAlternative {
                    reference: references[reference].label.clone(),
                    target: targets[target].label.clone(),
                    identity: pair_identity,
                });
            }
        }
    }
    alternatives.sort_by(|left, right| {
        left.reference
            .cmp(&right.reference)
            .then(right.identity.total_cmp(&left.identity))
            .then(left.target.cmp(&right.target))
    });
    alternatives
}

pub(super) fn optimal_pairs(scores: &[Vec<f64>], threshold: f64) -> Vec<(usize, usize)> {
    let rows = scores.len();
    let columns = scores.first().map_or(0, Vec::len);
    let size = rows.max(columns);
    if size == 0 {
        return Vec::new();
    }
    let costs = padded_costs(scores, threshold, size);
    hungarian(&costs)
        .into_iter()
        .enumerate()
        .filter(|(row, column)| {
            *row < rows && *column < columns && scores[*row][*column] >= threshold
        })
        .collect()
}

fn padded_costs(scores: &[Vec<f64>], threshold: f64, size: usize) -> Vec<Vec<f64>> {
    (0..size)
        .map(|row| {
            (0..size)
                .map(|column| {
                    let candidate = scores
                        .get(row)
                        .and_then(|values| values.get(column))
                        .copied()
                        .filter(|value| *value >= threshold);
                    let weight = match candidate {
                        Some(value) => value,
                        None => 0.0,
                    };
                    1.0 - weight
                })
                .collect()
        })
        .collect()
}

fn hungarian(costs: &[Vec<f64>]) -> Vec<usize> {
    let size = costs.len();
    let mut state = HungarianState::new(size);
    for row in 1..=size {
        state.augment(row, costs);
    }
    state.assignment()
}

struct HungarianState {
    row_potential: Vec<f64>,
    column_potential: Vec<f64>,
    column_row: Vec<usize>,
    predecessor: Vec<usize>,
}

impl HungarianState {
    fn new(size: usize) -> Self {
        Self {
            row_potential: vec![0.0; size + 1],
            column_potential: vec![0.0; size + 1],
            column_row: vec![0; size + 1],
            predecessor: vec![0; size + 1],
        }
    }

    fn augment(&mut self, row: usize, costs: &[Vec<f64>]) {
        self.column_row[0] = row;
        let mut frontier = SearchFrontier::new(costs.len());
        let mut column = 0usize;
        loop {
            frontier.used[column] = true;
            let current_row = self.column_row[column];
            let (delta, next) = self.update_frontier(current_row, column, costs, &mut frontier);
            self.apply_delta(delta, &mut frontier);
            column = next;
            if self.column_row[column] == 0 {
                break;
            }
        }
        self.install_path(column);
    }

    fn update_frontier(
        &mut self,
        row: usize,
        current_column: usize,
        costs: &[Vec<f64>],
        frontier: &mut SearchFrontier,
    ) -> (f64, usize) {
        let mut best = (f64::INFINITY, 0usize);
        for column in 1..=costs.len() {
            if frontier.used[column] {
                continue;
            }
            let reduced = costs[row - 1][column - 1]
                - self.row_potential[row]
                - self.column_potential[column];
            if reduced < frontier.minimum[column] {
                frontier.minimum[column] = reduced;
                self.predecessor[column] = current_column;
            }
            if frontier.minimum[column] < best.0 {
                best = (frontier.minimum[column], column);
            }
        }
        best
    }

    fn apply_delta(&mut self, delta: f64, frontier: &mut SearchFrontier) {
        for column in 0..self.column_row.len() {
            if frontier.used[column] {
                self.row_potential[self.column_row[column]] += delta;
                self.column_potential[column] -= delta;
            } else {
                frontier.minimum[column] -= delta;
            }
        }
    }

    fn install_path(&mut self, mut column: usize) {
        loop {
            let previous = self.predecessor[column];
            self.column_row[column] = self.column_row[previous];
            column = previous;
            if column == 0 {
                break;
            }
        }
    }

    fn assignment(&self) -> Vec<usize> {
        let mut result = vec![0; self.column_row.len() - 1];
        for column in 1..self.column_row.len() {
            result[self.column_row[column] - 1] = column - 1;
        }
        result
    }
}

struct SearchFrontier {
    minimum: Vec<f64>,
    used: Vec<bool>,
}

impl SearchFrontier {
    fn new(size: usize) -> Self {
        Self {
            minimum: vec![f64::INFINITY; size + 1],
            used: vec![false; size + 1],
        }
    }
}
