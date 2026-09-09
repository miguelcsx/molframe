//! Cache-sized dynamic-programming workspace for profile alignment.

use super::memory::alignment_peak_bytes;
use super::{GAP, Profile};
use crate::{MsaError, Scoring};

const NEGATIVE_INFINITY: i64 = i64::MIN / 4;
const MATCH_STATE: u8 = 0;
const LEFT_GAP_STATE: u8 = 1;
const RIGHT_GAP_STATE: u8 = 2;
const LEFT_EXTENDS: u8 = 1 << 2;
const RIGHT_EXTENDS: u8 = 1 << 3;

struct ColumnStats {
    offsets: Vec<u32>,
    symbols: Vec<u8>,
    counts: Vec<u32>,
    residues: Vec<u32>,
}

impl ColumnStats {
    fn new(profile: &Profile) -> Result<Self, MsaError> {
        let width = profile.width();
        let mut offsets = Vec::with_capacity(width.saturating_add(1));
        let mut symbols = Vec::new();
        let mut column_counts = Vec::new();
        let mut residues = Vec::with_capacity(width);
        let mut counts = [0_u32; 256];
        let mut touched = Vec::with_capacity(32);
        offsets.push(0);
        for column in 0..width {
            let mut residue_count = 0_u32;
            for (_, row) in &profile.rows {
                let symbol = row[column];
                if symbol == GAP {
                    continue;
                }
                let index = usize::from(symbol);
                if counts[index] == 0 {
                    touched.push(symbol);
                }
                counts[index] = counts[index]
                    .checked_add(1)
                    .ok_or(MsaError::DimensionOverflow)?;
                residue_count = residue_count
                    .checked_add(1)
                    .ok_or(MsaError::DimensionOverflow)?;
            }
            touched.sort_unstable();
            for symbol in touched.drain(..) {
                let index = usize::from(symbol);
                symbols.push(symbol);
                column_counts.push(counts[index]);
                counts[index] = 0;
            }
            offsets.push(u32::try_from(symbols.len()).map_err(|_| MsaError::DimensionOverflow)?);
            residues.push(residue_count);
        }
        Ok(Self {
            offsets,
            symbols,
            counts: column_counts,
            residues,
        })
    }

    fn range(&self, column: usize) -> std::ops::Range<usize> {
        self.offsets[column] as usize..self.offsets[column + 1] as usize
    }

    fn score(
        &self,
        left_column: usize,
        right: &Self,
        right_column: usize,
        scoring: Scoring,
    ) -> Result<i64, MsaError> {
        let total_pairs = i64::from(self.residues[left_column])
            .checked_mul(i64::from(right.residues[right_column]))
            .ok_or(MsaError::NumericOverflow)?;
        let mut matches = 0_i64;
        let left_range = self.range(left_column);
        let right_range = right.range(right_column);
        let mut left_index = left_range.start;
        let mut right_index = right_range.start;
        while left_index < left_range.end && right_index < right_range.end {
            match self.symbols[left_index].cmp(&right.symbols[right_index]) {
                std::cmp::Ordering::Less => left_index += 1,
                std::cmp::Ordering::Greater => right_index += 1,
                std::cmp::Ordering::Equal => {
                    matches = matches
                        .checked_add(
                            i64::from(self.counts[left_index])
                                .checked_mul(i64::from(right.counts[right_index]))
                                .ok_or(MsaError::NumericOverflow)?,
                        )
                        .ok_or(MsaError::NumericOverflow)?;
                    left_index += 1;
                    right_index += 1;
                }
            }
        }
        let mismatch_total = i64::from(scoring.mismatch_score)
            .checked_mul(total_pairs)
            .ok_or(MsaError::NumericOverflow)?;
        let match_bonus = i64::from(scoring.match_score)
            .checked_sub(i64::from(scoring.mismatch_score))
            .and_then(|bonus| bonus.checked_mul(matches))
            .ok_or(MsaError::NumericOverflow)?;
        mismatch_total
            .checked_add(match_bonus)
            .ok_or(MsaError::NumericOverflow)
    }
}

struct GapScores {
    open: Vec<i64>,
    extend: Vec<i64>,
}

impl GapScores {
    fn new(stats: &ColumnStats, opposite_rows: usize, scoring: Scoring) -> Result<Self, MsaError> {
        let opposite_rows =
            i64::try_from(opposite_rows).map_err(|_| MsaError::DimensionOverflow)?;
        let mut open = Vec::with_capacity(stats.residues.len());
        let mut extend = Vec::with_capacity(stats.residues.len());
        for &residues in &stats.residues {
            let pairs = i64::from(residues)
                .checked_mul(opposite_rows)
                .ok_or(MsaError::NumericOverflow)?;
            open.push(
                pairs
                    .checked_mul(i64::from(scoring.gap_open))
                    .ok_or(MsaError::NumericOverflow)?,
            );
            extend.push(
                pairs
                    .checked_mul(i64::from(scoring.gap_extend))
                    .ok_or(MsaError::NumericOverflow)?,
            );
        }
        Ok(Self { open, extend })
    }
}

struct ScoreRows {
    states: [Vec<i64>; 3],
}

impl ScoreRows {
    fn new(width: usize) -> Self {
        Self {
            states: std::array::from_fn(|_| vec![NEGATIVE_INFINITY; width]),
        }
    }

    fn get(&self, state: u8, column: usize) -> i64 {
        self.states[usize::from(state)][column]
    }

    fn set(&mut self, state: u8, column: usize, value: i64) {
        self.states[usize::from(state)][column] = value;
    }

    fn reset(&mut self) {
        for state in &mut self.states {
            state.fill(NEGATIVE_INFINITY);
        }
    }
}

struct PackedTrace {
    values: Vec<u8>,
}

impl PackedTrace {
    fn new(cells: usize) -> Result<Self, MsaError> {
        let bytes = Self::required_bytes(cells).ok_or(MsaError::DimensionOverflow)?;
        Ok(Self {
            values: vec![0; bytes],
        })
    }

    const fn required_bytes(cells: usize) -> Option<usize> {
        match cells.checked_add(1) {
            Some(value) => Some(value / 2),
            None => None,
        }
    }

    fn set(&mut self, index: usize, value: u8) {
        let byte = &mut self.values[index / 2];
        let shift = (index & 1) * 4;
        *byte = (*byte & !(0x0f << shift)) | ((value & 0x0f) << shift);
    }

    fn get(&self, index: usize) -> u8 {
        let shift = (index & 1) * 4;
        (self.values[index / 2] >> shift) & 0x0f
    }
}

struct CellContext<'a> {
    left_stats: &'a ColumnStats,
    right_stats: &'a ColumnStats,
    left_gaps: &'a GapScores,
    right_gaps: &'a GapScores,
    scoring: Scoring,
}

pub(in crate::msa) fn align(
    left: &Profile,
    right: &Profile,
    scoring: Scoring,
    memory_limit_bytes: usize,
) -> Result<Profile, MsaError> {
    let required = alignment_peak_bytes(left, right)?;
    if required > memory_limit_bytes {
        return Err(MsaError::MemoryLimit {
            required,
            limit: memory_limit_bytes,
        });
    }
    let left_stats = ColumnStats::new(left)?;
    let right_stats = ColumnStats::new(right)?;
    let left_gaps = GapScores::new(&left_stats, right.rows.len(), scoring)?;
    let right_gaps = GapScores::new(&right_stats, left.rows.len(), scoring)?;
    let trace_cells = left
        .width()
        .checked_mul(right.width())
        .ok_or(MsaError::DimensionOverflow)?;
    let mut trace = PackedTrace::new(trace_cells)?;
    let score_width = right
        .width()
        .checked_add(1)
        .ok_or(MsaError::DimensionOverflow)?;
    let mut previous = ScoreRows::new(score_width);
    let mut current = ScoreRows::new(score_width);
    let context = CellContext {
        left_stats: &left_stats,
        right_stats: &right_stats,
        left_gaps: &left_gaps,
        right_gaps: &right_gaps,
        scoring,
    };
    seed_horizontal(&mut previous, &right_gaps)?;
    let mut vertical = 0_i64;
    for row in 1..=left.width() {
        current.reset();
        let edge = if row == 1 {
            left_gaps.open[row - 1]
        } else {
            left_gaps.extend[row - 1]
        };
        vertical = add(vertical, edge)?;
        current.set(RIGHT_GAP_STATE, 0, vertical);
        for column in 1..=right.width() {
            context.fill(&previous, &mut current, row, column, &mut trace)?;
        }
        std::mem::swap(&mut previous, &mut current);
    }
    let (_, final_state) = best_state(&previous, right.width());
    traceback(left, right, &trace, final_state)
}

fn seed_horizontal(scores: &mut ScoreRows, gaps: &GapScores) -> Result<(), MsaError> {
    scores.set(MATCH_STATE, 0, 0);
    let mut score = 0_i64;
    for column in 1..=gaps.open.len() {
        let penalty = if column == 1 {
            gaps.open[column - 1]
        } else {
            gaps.extend[column - 1]
        };
        score = add(score, penalty)?;
        scores.set(LEFT_GAP_STATE, column, score);
    }
    Ok(())
}

impl CellContext<'_> {
    fn fill(
        &self,
        previous: &ScoreRows,
        current: &mut ScoreRows,
        row: usize,
        column: usize,
        trace: &mut PackedTrace,
    ) -> Result<(), MsaError> {
        let (diagonal, match_from) = best_state(previous, column - 1);
        current.set(
            MATCH_STATE,
            column,
            add(
                diagonal,
                self.left_stats
                    .score(row - 1, self.right_stats, column - 1, self.scoring)?,
            )?,
        );
        let (left_value, left_from) = prefer_first(
            add(
                current.get(MATCH_STATE, column - 1),
                self.right_gaps.open[column - 1],
            )?,
            MATCH_STATE,
            add(
                current.get(LEFT_GAP_STATE, column - 1),
                self.right_gaps.extend[column - 1],
            )?,
            LEFT_GAP_STATE,
        );
        current.set(LEFT_GAP_STATE, column, left_value);
        let (right_value, right_from) = prefer_first(
            add(
                previous.get(MATCH_STATE, column),
                self.left_gaps.open[row - 1],
            )?,
            MATCH_STATE,
            add(
                previous.get(RIGHT_GAP_STATE, column),
                self.left_gaps.extend[row - 1],
            )?,
            RIGHT_GAP_STATE,
        );
        current.set(RIGHT_GAP_STATE, column, right_value);
        let value = match_from
            | if left_from == LEFT_GAP_STATE {
                LEFT_EXTENDS
            } else {
                0
            }
            | if right_from == RIGHT_GAP_STATE {
                RIGHT_EXTENDS
            } else {
                0
            };
        trace.set(
            trace_index(row, column, self.right_stats.residues.len()),
            value,
        );
        Ok(())
    }
}

fn trace_index(row: usize, column: usize, columns: usize) -> usize {
    (row - 1) * columns + column - 1
}

fn traceback(
    left: &Profile,
    right: &Profile,
    trace: &PackedTrace,
    mut state: u8,
) -> Result<Profile, MsaError> {
    let mut row = left.width();
    let mut column = right.width();
    let capacity = row.checked_add(column).ok_or(MsaError::DimensionOverflow)?;
    let mut left_columns = Vec::with_capacity(capacity);
    let mut right_columns = Vec::with_capacity(capacity);
    while row > 0 || column > 0 {
        match state {
            MATCH_STATE if row > 0 && column > 0 => {
                let value = trace.get(trace_index(row, column, right.width()));
                row -= 1;
                column -= 1;
                left_columns.push(Some(row));
                right_columns.push(Some(column));
                state = value & 0b11;
            }
            LEFT_GAP_STATE if column > 0 => {
                let previous = if row == 0 {
                    if column == 1 {
                        MATCH_STATE
                    } else {
                        LEFT_GAP_STATE
                    }
                } else {
                    let value = trace.get(trace_index(row, column, right.width()));
                    if value & LEFT_EXTENDS == 0 {
                        MATCH_STATE
                    } else {
                        LEFT_GAP_STATE
                    }
                };
                column -= 1;
                left_columns.push(None);
                right_columns.push(Some(column));
                state = previous;
            }
            RIGHT_GAP_STATE if row > 0 => {
                let previous = if column == 0 {
                    if row == 1 {
                        MATCH_STATE
                    } else {
                        RIGHT_GAP_STATE
                    }
                } else {
                    let value = trace.get(trace_index(row, column, right.width()));
                    if value & RIGHT_EXTENDS == 0 {
                        MATCH_STATE
                    } else {
                        RIGHT_GAP_STATE
                    }
                };
                row -= 1;
                left_columns.push(Some(row));
                right_columns.push(None);
                state = previous;
            }
            _ => return Err(MsaError::NumericOverflow),
        }
    }
    left_columns.reverse();
    right_columns.reverse();
    let mut rows = project(left, &left_columns);
    rows.extend(project(right, &right_columns));
    Ok(Profile { rows })
}

fn project(profile: &Profile, columns: &[Option<usize>]) -> Vec<(usize, Vec<u8>)> {
    profile
        .rows
        .iter()
        .map(|(index, row)| {
            let projected = columns
                .iter()
                .map(|source| match source {
                    Some(source) => row[*source],
                    None => GAP,
                })
                .collect();
            (*index, projected)
        })
        .collect()
}

fn best_state(scores: &ScoreRows, column: usize) -> (i64, u8) {
    let mut best = (scores.get(MATCH_STATE, column), MATCH_STATE);
    best = prefer_first(
        best.0,
        best.1,
        scores.get(RIGHT_GAP_STATE, column),
        RIGHT_GAP_STATE,
    );
    prefer_first(
        best.0,
        best.1,
        scores.get(LEFT_GAP_STATE, column),
        LEFT_GAP_STATE,
    )
}

fn prefer_first(first: i64, first_state: u8, second: i64, second_state: u8) -> (i64, u8) {
    if first >= second {
        (first, first_state)
    } else {
        (second, second_state)
    }
}

fn add(value: i64, addition: i64) -> Result<i64, MsaError> {
    if value <= NEGATIVE_INFINITY {
        Ok(NEGATIVE_INFINITY)
    } else {
        value.checked_add(addition).ok_or(MsaError::NumericOverflow)
    }
}

#[cfg(test)]
#[path = "workspace_tests.rs"]
mod tests;
