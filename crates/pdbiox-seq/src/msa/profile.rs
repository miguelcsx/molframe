//! Affine profile-profile alignment and sum-of-pairs scoring.

use super::MsaError;
use crate::{Score, Scoring};

pub(super) const GAP: u8 = b'-';
const NEGATIVE_INFINITY: i64 = i64::MIN / 4;
const MATCH_STATE: u8 = 0;
const LEFT_GAP_STATE: u8 = 1;
const RIGHT_GAP_STATE: u8 = 2;

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct Profile {
    pub(super) rows: Vec<(usize, Vec<u8>)>,
}

impl Profile {
    pub(super) fn singleton(index: usize, sequence: &[u8]) -> Self {
        Self {
            rows: vec![(index, sequence.to_vec())],
        }
    }

    pub(super) fn width(&self) -> usize {
        match self.rows.first() {
            Some((_, row)) => row.len(),
            None => 0,
        }
    }

    fn column(&self, index: usize) -> impl Iterator<Item = u8> + '_ {
        self.rows.iter().map(move |(_, row)| row[index])
    }

    pub(super) fn without_row(&self, row_index: usize) -> Self {
        let mut rows = self
            .rows
            .iter()
            .filter(|(index, _)| *index != row_index)
            .cloned()
            .collect::<Vec<_>>();
        remove_all_gap_columns(&mut rows);
        Self { rows }
    }

    pub(super) fn score(&self, scoring: Scoring) -> Result<i64, MsaError> {
        let mut score = 0_i64;
        for left in 0..self.rows.len() {
            for right in left + 1..self.rows.len() {
                score = score
                    .checked_add(pair_score(
                        &self.rows[left].1,
                        &self.rows[right].1,
                        scoring,
                    )?)
                    .ok_or(MsaError::NumericOverflow)?;
            }
        }
        Ok(score)
    }
}

struct Tables {
    width: usize,
    values: [Vec<i64>; 3],
    trace: [Vec<u8>; 3],
}

impl Tables {
    fn index(&self, row: usize, column: usize) -> usize {
        row * self.width + column
    }
}

pub(super) fn align(
    left: &Profile,
    right: &Profile,
    scoring: Scoring,
) -> Result<Profile, MsaError> {
    let rows = left
        .width()
        .checked_add(1)
        .ok_or(MsaError::DimensionOverflow)?;
    let width = right
        .width()
        .checked_add(1)
        .ok_or(MsaError::DimensionOverflow)?;
    let size = rows.checked_mul(width).ok_or(MsaError::DimensionOverflow)?;
    let mut tables = Tables {
        width,
        values: std::array::from_fn(|_| vec![NEGATIVE_INFINITY; size]),
        trace: std::array::from_fn(|_| vec![MATCH_STATE; size]),
    };
    tables.values[usize::from(MATCH_STATE)][0] = 0;
    seed_edges(&mut tables, left, right, scoring)?;
    for row in 1..rows {
        for column in 1..width {
            fill(&mut tables, left, right, scoring, row, column)?;
        }
    }
    traceback(&tables, left, right)
}

fn seed_edges(
    tables: &mut Tables,
    left: &Profile,
    right: &Profile,
    scoring: Scoring,
) -> Result<(), MsaError> {
    let mut score = 0_i64;
    for row in 1..=left.width() {
        score = score
            .checked_add(gap_column_score(
                left,
                row - 1,
                right.rows.len(),
                row == 1,
                scoring,
            )?)
            .ok_or(MsaError::NumericOverflow)?;
        let index = tables.index(row, 0);
        tables.values[usize::from(RIGHT_GAP_STATE)][index] = score;
        tables.trace[usize::from(RIGHT_GAP_STATE)][index] = if row == 1 {
            MATCH_STATE
        } else {
            RIGHT_GAP_STATE
        };
    }
    score = 0;
    for column in 1..=right.width() {
        score = score
            .checked_add(gap_column_score(
                right,
                column - 1,
                left.rows.len(),
                column == 1,
                scoring,
            )?)
            .ok_or(MsaError::NumericOverflow)?;
        let index = tables.index(0, column);
        tables.values[usize::from(LEFT_GAP_STATE)][index] = score;
        tables.trace[usize::from(LEFT_GAP_STATE)][index] = if column == 1 {
            MATCH_STATE
        } else {
            LEFT_GAP_STATE
        };
    }
    Ok(())
}

fn fill(
    tables: &mut Tables,
    left: &Profile,
    right: &Profile,
    scoring: Scoring,
    row: usize,
    column: usize,
) -> Result<(), MsaError> {
    let here = tables.index(row, column);
    let diagonal = tables.index(row - 1, column - 1);
    let (value, state) = best_state(&tables.values, diagonal);
    tables.values[usize::from(MATCH_STATE)][here] = add(
        value,
        column_score(left, row - 1, right, column - 1, scoring)?,
    )?;
    tables.trace[usize::from(MATCH_STATE)][here] = state;

    let previous = tables.index(row, column - 1);
    let open = add(
        tables.values[usize::from(MATCH_STATE)][previous],
        gap_column_score(right, column - 1, left.rows.len(), true, scoring)?,
    )?;
    let extend = add(
        tables.values[usize::from(LEFT_GAP_STATE)][previous],
        gap_column_score(right, column - 1, left.rows.len(), false, scoring)?,
    )?;
    let (value, state) = prefer_first(open, MATCH_STATE, extend, LEFT_GAP_STATE);
    tables.values[usize::from(LEFT_GAP_STATE)][here] = value;
    tables.trace[usize::from(LEFT_GAP_STATE)][here] = state;

    let previous = tables.index(row - 1, column);
    let open = add(
        tables.values[usize::from(MATCH_STATE)][previous],
        gap_column_score(left, row - 1, right.rows.len(), true, scoring)?,
    )?;
    let extend = add(
        tables.values[usize::from(RIGHT_GAP_STATE)][previous],
        gap_column_score(left, row - 1, right.rows.len(), false, scoring)?,
    )?;
    let (value, state) = prefer_first(open, MATCH_STATE, extend, RIGHT_GAP_STATE);
    tables.values[usize::from(RIGHT_GAP_STATE)][here] = value;
    tables.trace[usize::from(RIGHT_GAP_STATE)][here] = state;
    Ok(())
}

fn traceback(tables: &Tables, left: &Profile, right: &Profile) -> Result<Profile, MsaError> {
    let mut row = left.width();
    let mut column = right.width();
    let final_index = tables.index(row, column);
    let (_, mut state) = best_state(&tables.values, final_index);
    let capacity = row.checked_add(column).ok_or(MsaError::DimensionOverflow)?;
    let mut left_columns = Vec::with_capacity(capacity);
    let mut right_columns = Vec::with_capacity(capacity);
    while row > 0 || column > 0 {
        let here = tables.index(row, column);
        let previous = tables.trace[usize::from(state)][here];
        match state {
            MATCH_STATE => {
                row -= 1;
                column -= 1;
                left_columns.push(Some(row));
                right_columns.push(Some(column));
            }
            LEFT_GAP_STATE => {
                column -= 1;
                left_columns.push(None);
                right_columns.push(Some(column));
            }
            RIGHT_GAP_STATE => {
                row -= 1;
                left_columns.push(Some(row));
                right_columns.push(None);
            }
            _ => break,
        }
        state = previous;
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

fn best_state(values: &[Vec<i64>; 3], index: usize) -> (i64, u8) {
    let mut best = (values[usize::from(MATCH_STATE)][index], MATCH_STATE);
    best = prefer_first(
        best.0,
        best.1,
        values[usize::from(RIGHT_GAP_STATE)][index],
        RIGHT_GAP_STATE,
    );
    prefer_first(
        best.0,
        best.1,
        values[usize::from(LEFT_GAP_STATE)][index],
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

fn column_score(
    left: &Profile,
    left_column: usize,
    right: &Profile,
    right_column: usize,
    scoring: Scoring,
) -> Result<i64, MsaError> {
    let mut total = 0_i64;
    for first in left.column(left_column) {
        for second in right.column(right_column) {
            total = total
                .checked_add(match (first, second) {
                    (GAP, _) | (_, GAP) => 0,
                    _ => i64::from(scoring.score(first, second)),
                })
                .ok_or(MsaError::NumericOverflow)?;
        }
    }
    Ok(total)
}

fn gap_column_score(
    profile: &Profile,
    column: usize,
    opposite_rows: usize,
    opening: bool,
    scoring: Scoring,
) -> Result<i64, MsaError> {
    let residues = profile
        .column(column)
        .filter(|symbol| *symbol != GAP)
        .count();
    let penalty = if opening {
        scoring.gap_open
    } else {
        scoring.gap_extend
    };
    let residues = i64::try_from(residues).map_err(|_| MsaError::DimensionOverflow)?;
    let opposite_rows = i64::try_from(opposite_rows).map_err(|_| MsaError::DimensionOverflow)?;
    i64::from(penalty)
        .checked_mul(residues)
        .and_then(|value| value.checked_mul(opposite_rows))
        .ok_or(MsaError::NumericOverflow)
}

fn pair_score(left: &[u8], right: &[u8], scoring: Scoring) -> Result<i64, MsaError> {
    let mut total = 0_i64;
    let mut gap_state = None;
    for (&first, &second) in left.iter().zip(right) {
        match (first == GAP, second == GAP) {
            (false, false) => {
                total = total
                    .checked_add(i64::from(scoring.score(first, second)))
                    .ok_or(MsaError::NumericOverflow)?;
                gap_state = None;
            }
            (true, true) => {}
            (left_gap, _) => {
                let current = Some(left_gap);
                let penalty = if gap_state == current {
                    scoring.gap_extend
                } else {
                    scoring.gap_open
                };
                total = total
                    .checked_add(i64::from(penalty))
                    .ok_or(MsaError::NumericOverflow)?;
                gap_state = current;
            }
        }
    }
    Ok(total)
}

fn remove_all_gap_columns(rows: &mut [(usize, Vec<u8>)]) {
    let width = match rows.first() {
        Some((_, row)) => row.len(),
        None => return,
    };
    let keep = (0..width)
        .filter(|column| rows.iter().any(|(_, row)| row[*column] != GAP))
        .collect::<Vec<_>>();
    for (_, row) in rows {
        *row = keep.iter().map(|column| row[*column]).collect();
    }
}
