//! Linear-space affine global traceback by divide and conquer.

use super::small::small_trace;
use super::types::{AlignError, Alignment, Column};
use crate::Score;

pub(super) const NEGATIVE_INFINITY: i64 = i64::MIN / 4;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum State {
    Match,
    LeftGap,
    RightGap,
}

impl State {
    pub(super) const ALL: [Self; 3] = [Self::Match, Self::RightGap, Self::LeftGap];

    pub(super) const fn index(self) -> usize {
        match self {
            Self::Match => 0,
            Self::LeftGap => 1,
            Self::RightGap => 2,
        }
    }
}

#[derive(Clone)]
struct RowScores {
    states: [Vec<i64>; 3],
}

pub(super) struct LinearProblem<'a, S> {
    pub(super) scorer: &'a S,
    pub(super) gap_open: i32,
    pub(super) gap_extend: i32,
}

#[derive(Clone, Copy)]
pub(super) struct TraceBoundary {
    pub(super) start: State,
    pub(super) end: State,
    pub(super) row_offset: usize,
    pub(super) column_offset: usize,
}

impl RowScores {
    fn new(width: usize) -> Self {
        Self {
            states: std::array::from_fn(|_| vec![NEGATIVE_INFINITY; width]),
        }
    }

    fn get(&self, state: State, column: usize) -> i64 {
        self.states[state.index()][column]
    }

    fn set(&mut self, state: State, column: usize, value: i64) {
        self.states[state.index()][column] = value;
    }
}

pub(super) fn global_linear<S: Score>(
    left: &[u8],
    right: &[u8],
    scorer: &S,
    gap_open: i32,
    gap_extend: i32,
) -> Result<Alignment, AlignError> {
    if right.len() <= left.len() {
        solve(left, right, scorer, gap_open, gap_extend, false)
    } else {
        solve(
            right,
            left,
            &ReversedScore(scorer),
            gap_open,
            gap_extend,
            true,
        )
    }
}

fn solve<S: Score>(
    rows: &[u8],
    columns: &[u8],
    scorer: &S,
    gap_open: i32,
    gap_extend: i32,
    transpose: bool,
) -> Result<Alignment, AlignError> {
    let final_row = forward(rows, columns, scorer, gap_open, gap_extend, State::Match)?;
    let (score, end) = best_end(&final_row, columns.len());
    let problem = LinearProblem {
        scorer,
        gap_open,
        gap_extend,
    };
    let columns_result = divide(
        rows,
        columns,
        &problem,
        TraceBoundary {
            start: State::Match,
            end,
            row_offset: 0,
            column_offset: 0,
        },
    )?;
    let Some(columns) = columns_result else {
        return Err(AlignError::NoAlignmentPath);
    };
    let mut alignment = Alignment {
        score: i32::try_from(score).map_err(|_| AlignError::NumericOverflow)?,
        columns,
    };
    if transpose {
        for column in &mut alignment.columns {
            std::mem::swap(&mut column.left, &mut column.right);
        }
    }
    Ok(alignment)
}

#[derive(Clone, Copy)]
struct ReversedScore<'a, S>(&'a S);

impl<S: Score> Score for ReversedScore<'_, S> {
    fn score(&self, first: u8, second: u8) -> i32 {
        self.0.score(second, first)
    }
}

fn divide<S: Score>(
    rows: &[u8],
    columns: &[u8],
    problem: &LinearProblem<'_, S>,
    boundary: TraceBoundary,
) -> Result<Option<Vec<Column>>, AlignError> {
    if rows.len() <= 1 || columns.len() <= 1 {
        return small_trace(rows, columns, problem, boundary);
    }
    let middle = rows.len() / 2;
    let forward_scores = forward(
        &rows[..middle],
        columns,
        problem.scorer,
        problem.gap_open,
        problem.gap_extend,
        boundary.start,
    )?;
    let backward_scores = backward(
        &rows[middle..],
        columns,
        problem.scorer,
        problem.gap_open,
        problem.gap_extend,
        boundary.end,
    )?;
    let Some((split_column, split_state)) = split(&forward_scores, &backward_scores)? else {
        return Ok(None);
    };
    let Some(mut first) = divide(
        &rows[..middle],
        &columns[..split_column],
        problem,
        TraceBoundary {
            end: split_state,
            ..boundary
        },
    )?
    else {
        return Ok(None);
    };
    let row_offset = boundary
        .row_offset
        .checked_add(middle)
        .ok_or(AlignError::NumericOverflow)?;
    let column_offset = boundary
        .column_offset
        .checked_add(split_column)
        .ok_or(AlignError::NumericOverflow)?;
    let Some(second) = divide(
        &rows[middle..],
        &columns[split_column..],
        problem,
        TraceBoundary {
            start: split_state,
            row_offset,
            column_offset,
            ..boundary
        },
    )?
    else {
        return Ok(None);
    };
    first.extend(second);
    Ok(Some(first))
}

fn forward<S: Score>(
    rows: &[u8],
    columns: &[u8],
    scorer: &S,
    gap_open: i32,
    gap_extend: i32,
    start: State,
) -> Result<RowScores, AlignError> {
    let width = columns
        .len()
        .checked_add(1)
        .ok_or(AlignError::NumericOverflow)?;
    let mut previous = RowScores::new(width);
    previous.set(start, 0, 0);
    for column in 1..width {
        let value = best_two(
            add(previous.get(State::Match, column - 1), gap_open)?,
            add(previous.get(State::LeftGap, column - 1), gap_extend)?,
        );
        previous.set(State::LeftGap, column, value);
    }
    for row_symbol in rows.iter().copied() {
        let mut current = RowScores::new(width);
        current.set(
            State::RightGap,
            0,
            best_two(
                add(previous.get(State::Match, 0), gap_open)?,
                add(previous.get(State::RightGap, 0), gap_extend)?,
            ),
        );
        for (column_index, column_symbol) in columns.iter().copied().enumerate() {
            let column = column_index + 1;
            let diagonal = best_three(
                previous.get(State::Match, column - 1),
                previous.get(State::RightGap, column - 1),
                previous.get(State::LeftGap, column - 1),
            );
            current.set(
                State::Match,
                column,
                add(diagonal, scorer.score(row_symbol, column_symbol))?,
            );
            current.set(
                State::RightGap,
                column,
                best_two(
                    add(previous.get(State::Match, column), gap_open)?,
                    add(previous.get(State::RightGap, column), gap_extend)?,
                ),
            );
            current.set(
                State::LeftGap,
                column,
                best_two(
                    add(current.get(State::Match, column - 1), gap_open)?,
                    add(current.get(State::LeftGap, column - 1), gap_extend)?,
                ),
            );
        }
        previous = current;
    }
    Ok(previous)
}

fn backward<S: Score>(
    rows: &[u8],
    columns: &[u8],
    scorer: &S,
    gap_open: i32,
    gap_extend: i32,
    end: State,
) -> Result<RowScores, AlignError> {
    let width = columns
        .len()
        .checked_add(1)
        .ok_or(AlignError::NumericOverflow)?;
    let mut next = RowScores::new(width);
    next.set(end, columns.len(), 0);
    for column in (0..columns.len()).rev() {
        next.set(
            State::Match,
            column,
            add(next.get(State::LeftGap, column + 1), gap_open)?,
        );
        next.set(
            State::LeftGap,
            column,
            add(next.get(State::LeftGap, column + 1), gap_extend)?,
        );
    }
    for row in (0..rows.len()).rev() {
        let mut current = RowScores::new(width);
        current.set(
            State::Match,
            columns.len(),
            add(next.get(State::RightGap, columns.len()), gap_open)?,
        );
        current.set(
            State::RightGap,
            columns.len(),
            add(next.get(State::RightGap, columns.len()), gap_extend)?,
        );
        for column in (0..columns.len()).rev() {
            let diagonal = add(
                next.get(State::Match, column + 1),
                scorer.score(rows[row], columns[column]),
            )?;
            current.set(
                State::Match,
                column,
                best_three(
                    diagonal,
                    add(next.get(State::RightGap, column), gap_open)?,
                    add(current.get(State::LeftGap, column + 1), gap_open)?,
                ),
            );
            current.set(
                State::RightGap,
                column,
                best_two(
                    diagonal,
                    add(next.get(State::RightGap, column), gap_extend)?,
                ),
            );
            current.set(
                State::LeftGap,
                column,
                best_two(
                    diagonal,
                    add(current.get(State::LeftGap, column + 1), gap_extend)?,
                ),
            );
        }
        next = current;
    }
    Ok(next)
}

fn split(forward: &RowScores, backward: &RowScores) -> Result<Option<(usize, State)>, AlignError> {
    let width = forward.states[0].len();
    let mut best: Option<(i64, usize, State)> = None;
    for column in 0..width {
        for state in State::ALL {
            let score = add_exact(forward.get(state, column), backward.get(state, column))?;
            match best {
                Some((current, _, _)) if score <= current => {}
                _ => best = Some((score, column, state)),
            }
        }
    }
    Ok(best
        .filter(|(score, _, _)| *score > NEGATIVE_INFINITY)
        .map(|(_, column, state)| (column, state)))
}

fn best_end(row: &RowScores, column: usize) -> (i64, State) {
    let mut best = (row.get(State::Match, column), State::Match);
    for state in [State::RightGap, State::LeftGap] {
        let candidate = row.get(state, column);
        if candidate > best.0 {
            best = (candidate, state);
        }
    }
    best
}

pub(super) fn add(value: i64, addition: i32) -> Result<i64, AlignError> {
    if value <= NEGATIVE_INFINITY {
        Ok(NEGATIVE_INFINITY)
    } else {
        value
            .checked_add(i64::from(addition))
            .ok_or(AlignError::NumericOverflow)
    }
}

fn add_exact(first: i64, second: i64) -> Result<i64, AlignError> {
    if first <= NEGATIVE_INFINITY || second <= NEGATIVE_INFINITY {
        Ok(NEGATIVE_INFINITY)
    } else {
        first.checked_add(second).ok_or(AlignError::NumericOverflow)
    }
}

const fn best_two(first: i64, second: i64) -> i64 {
    if first >= second { first } else { second }
}

const fn best_three(first: i64, second: i64, third: i64) -> i64 {
    best_two(best_two(first, second), third)
}

#[cfg(test)]
#[path = "linear_tests.rs"]
mod tests;
