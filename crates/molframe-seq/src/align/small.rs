//! Bounded affine traceback base case for the linear-space solver.

use super::linear::{LinearProblem, NEGATIVE_INFINITY, State, TraceBoundary, add};
use super::types::{AlignError, Column};
use crate::Score;

struct Workspace {
    width: usize,
    scores: [Vec<i64>; 3],
    trace: [Vec<Option<State>>; 3],
}

pub(super) fn small_trace<S: Score>(
    rows: &[u8],
    columns: &[u8],
    problem: &LinearProblem<'_, S>,
    boundary: TraceBoundary,
) -> Result<Option<Vec<Column>>, AlignError> {
    let width = columns
        .len()
        .checked_add(1)
        .ok_or(AlignError::NumericOverflow)?;
    let height = rows
        .len()
        .checked_add(1)
        .ok_or(AlignError::NumericOverflow)?;
    let size = height
        .checked_mul(width)
        .ok_or(AlignError::NumericOverflow)?;
    let mut workspace = Workspace {
        width,
        scores: std::array::from_fn(|_| vec![NEGATIVE_INFINITY; size]),
        trace: std::array::from_fn(|_| vec![None; size]),
    };
    workspace.scores[boundary.start.index()][0] = 0;
    fill(rows, columns, problem, &mut workspace)?;
    traceback(rows.len(), columns.len(), boundary, &workspace)
}

fn fill<S: Score>(
    rows: &[u8],
    columns: &[u8],
    problem: &LinearProblem<'_, S>,
    workspace: &mut Workspace,
) -> Result<(), AlignError> {
    for row in 0..=rows.len() {
        for column in 0..=columns.len() {
            if row == 0 && column == 0 {
                continue;
            }
            fill_cell(row, column, rows, columns, problem, workspace)?;
        }
    }
    Ok(())
}

fn fill_cell<S: Score>(
    row: usize,
    column: usize,
    rows: &[u8],
    columns: &[u8],
    problem: &LinearProblem<'_, S>,
    workspace: &mut Workspace,
) -> Result<(), AlignError> {
    let index = row * workspace.width + column;
    if row > 0 && column > 0 {
        let previous = (row - 1) * workspace.width + column - 1;
        let (value, state) = best_predecessor(&workspace.scores, previous);
        workspace.scores[State::Match.index()][index] = add(
            value,
            problem.scorer.score(rows[row - 1], columns[column - 1]),
        )?;
        workspace.trace[State::Match.index()][index] = Some(state);
    }
    if row > 0 {
        let previous = (row - 1) * workspace.width + column;
        let choices = [
            (
                add(
                    workspace.scores[State::Match.index()][previous],
                    problem.gap_open,
                )?,
                State::Match,
            ),
            (
                add(
                    workspace.scores[State::RightGap.index()][previous],
                    problem.gap_extend,
                )?,
                State::RightGap,
            ),
        ];
        let (value, state) = best_choice(&choices);
        workspace.scores[State::RightGap.index()][index] = value;
        workspace.trace[State::RightGap.index()][index] = Some(state);
    }
    if column > 0 {
        let previous = row * workspace.width + column - 1;
        let choices = [
            (
                add(
                    workspace.scores[State::Match.index()][previous],
                    problem.gap_open,
                )?,
                State::Match,
            ),
            (
                add(
                    workspace.scores[State::LeftGap.index()][previous],
                    problem.gap_extend,
                )?,
                State::LeftGap,
            ),
        ];
        let (value, state) = best_choice(&choices);
        workspace.scores[State::LeftGap.index()][index] = value;
        workspace.trace[State::LeftGap.index()][index] = Some(state);
    }
    Ok(())
}

fn traceback(
    mut row: usize,
    mut column: usize,
    boundary: TraceBoundary,
    workspace: &Workspace,
) -> Result<Option<Vec<Column>>, AlignError> {
    let mut state = boundary.end;
    if workspace.scores[state.index()][row * workspace.width + column] <= NEGATIVE_INFINITY {
        return Ok(None);
    }
    let capacity = row.checked_add(column).ok_or(AlignError::NumericOverflow)?;
    let mut result = Vec::with_capacity(capacity);
    while row > 0 || column > 0 {
        let Some(previous) = workspace.trace[state.index()][row * workspace.width + column] else {
            return Ok(None);
        };
        result.push(previous_column(&mut row, &mut column, state, boundary));
        state = previous;
    }
    result.reverse();
    Ok(Some(result))
}

fn previous_column(
    row: &mut usize,
    column: &mut usize,
    state: State,
    boundary: TraceBoundary,
) -> Column {
    match state {
        State::Match => {
            *row -= 1;
            *column -= 1;
            Column {
                left: Some(boundary.row_offset + *row),
                right: Some(boundary.column_offset + *column),
            }
        }
        State::RightGap => {
            *row -= 1;
            Column {
                left: Some(boundary.row_offset + *row),
                right: None,
            }
        }
        State::LeftGap => {
            *column -= 1;
            Column {
                left: None,
                right: Some(boundary.column_offset + *column),
            }
        }
    }
}

fn best_predecessor(scores: &[Vec<i64>; 3], index: usize) -> (i64, State) {
    let choices = State::ALL.map(|state| (scores[state.index()][index], state));
    best_choice(&choices)
}

fn best_choice(choices: &[(i64, State)]) -> (i64, State) {
    let mut best = choices[0];
    for choice in &choices[1..] {
        if choice.0 > best.0 {
            best = *choice;
        }
    }
    best
}
