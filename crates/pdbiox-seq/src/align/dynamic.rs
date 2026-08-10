//! Shared Gotoh matrices and deterministic traceback.

use super::types::{AlignError, Alignment, Column};
use crate::scoring::Score;

const NEG: i64 = i64::MIN / 4;
pub(super) const FROM_M: u8 = 0;
pub(super) const FROM_IX: u8 = 1;
pub(super) const FROM_IY: u8 = 2;
pub(super) const STOP: u8 = 3;

#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) enum Mode {
    Global,
    Local,
    SemiGlobal,
}

struct Tables {
    width: usize,
    m: Vec<i64>,
    ix: Vec<i64>,
    iy: Vec<i64>,
    from_m: Vec<u8>,
    from_ix: Vec<u8>,
    from_iy: Vec<u8>,
}

impl Tables {
    fn at(&self, row: usize, col: usize) -> usize {
        row * self.width + col
    }
}

pub(super) fn run<S: Score>(
    left: &[u8],
    right: &[u8],
    scorer: &S,
    gap_open: i32,
    gap_extend: i32,
    mode: Mode,
    band: Option<usize>,
) -> Result<Alignment, AlignError> {
    let rows = left
        .len()
        .checked_add(1)
        .ok_or(AlignError::NumericOverflow)?;
    let width = right
        .len()
        .checked_add(1)
        .ok_or(AlignError::NumericOverflow)?;
    let mut tables = seed(rows, width, gap_open, gap_extend, mode)?;
    let cell = CellParameters {
        left,
        right,
        scorer,
        gap_open,
        gap_extend,
        mode,
    };

    for i in 1..rows {
        for j in 1..width {
            if let Some(band) = band
                && i.abs_diff(j) > band
            {
                continue;
            }
            fill_cell(&mut tables, &cell, i, j)?;
        }
    }

    let (start_row, start_col, start_state, score) = start(&tables, mode);
    if score <= NEG {
        return Err(AlignError::NoAlignmentPath);
    }
    let columns = trace(&tables, mode, start_row, start_col, start_state);
    Ok(Alignment {
        score: i32::try_from(score).map_err(|_| AlignError::NumericOverflow)?,
        columns,
    })
}

fn seed(
    rows: usize,
    width: usize,
    gap_open: i32,
    gap_extend: i32,
    mode: Mode,
) -> Result<Tables, AlignError> {
    let size = rows.checked_mul(width).ok_or(AlignError::NumericOverflow)?;
    let mut tables = Tables {
        width,
        m: vec![NEG; size],
        ix: vec![NEG; size],
        iy: vec![NEG; size],
        from_m: vec![STOP; size],
        from_ix: vec![STOP; size],
        from_iy: vec![STOP; size],
    };
    tables.m[0] = 0;
    match mode {
        Mode::Global => {
            let mut penalty = i64::from(gap_open);
            for i in 1..rows {
                let index = tables.at(i, 0);
                tables.ix[index] = penalty;
                if i + 1 < rows {
                    penalty = penalty
                        .checked_add(i64::from(gap_extend))
                        .ok_or(AlignError::NumericOverflow)?;
                }
            }
            let mut penalty = i64::from(gap_open);
            for j in 1..width {
                let index = tables.at(0, j);
                tables.iy[index] = penalty;
                if j + 1 < width {
                    penalty = penalty
                        .checked_add(i64::from(gap_extend))
                        .ok_or(AlignError::NumericOverflow)?;
                }
            }
        }
        Mode::Local | Mode::SemiGlobal => {
            for i in 0..rows {
                let index = tables.at(i, 0);
                tables.m[index] = 0;
                tables.ix[index] = 0;
            }
            for j in 0..width {
                let index = tables.at(0, j);
                tables.m[index] = 0;
                tables.iy[index] = 0;
            }
        }
    }
    Ok(tables)
}

struct CellParameters<'a, S> {
    left: &'a [u8],
    right: &'a [u8],
    scorer: &'a S,
    gap_open: i32,
    gap_extend: i32,
    mode: Mode,
}

fn fill_cell<S: Score>(
    tables: &mut Tables,
    parameters: &CellParameters<'_, S>,
    i: usize,
    j: usize,
) -> Result<(), AlignError> {
    let here = tables.at(i, j);
    let diagonal = tables.at(i - 1, j - 1);
    let (diagonal_best, diagonal_from) =
        best3(tables.m[diagonal], tables.ix[diagonal], tables.iy[diagonal]);
    let substitution = parameters
        .scorer
        .score(parameters.left[i - 1], parameters.right[j - 1]);
    let mut match_value = add(diagonal_best, substitution)?;
    let mut match_from = diagonal_from;
    if parameters.mode == Mode::Local && match_value < 0 {
        match_value = 0;
        match_from = STOP;
    }
    tables.m[here] = match_value;
    tables.from_m[here] = match_from;

    let above = tables.at(i - 1, j);
    let (vertical_gap_value, vertical_gap_from) = better(
        add(tables.m[above], parameters.gap_open)?,
        FROM_M,
        add(tables.ix[above], parameters.gap_extend)?,
        FROM_IX,
    );
    tables.ix[here] = vertical_gap_value;
    tables.from_ix[here] = vertical_gap_from;

    let leftward = tables.at(i, j - 1);
    let (horizontal_gap_value, horizontal_gap_from) = better(
        add(tables.m[leftward], parameters.gap_open)?,
        FROM_M,
        add(tables.iy[leftward], parameters.gap_extend)?,
        FROM_IY,
    );
    tables.iy[here] = horizontal_gap_value;
    tables.from_iy[here] = horizontal_gap_from;
    Ok(())
}

fn add(value: i64, addition: i32) -> Result<i64, AlignError> {
    if value <= NEG {
        Ok(NEG)
    } else {
        value
            .checked_add(i64::from(addition))
            .ok_or(AlignError::NumericOverflow)
    }
}

/// Picks the best of the three diagonal predecessors, preferring M, then Ix.
fn best3(m: i64, ix: i64, iy: i64) -> (i64, u8) {
    let (value, from) = better(m, FROM_M, ix, FROM_IX);
    better(value, from, iy, FROM_IY)
}

/// Returns the higher-scoring of two options, keeping the first on a tie.
fn better(first: i64, first_from: u8, second: i64, second_from: u8) -> (i64, u8) {
    if first >= second {
        (first, first_from)
    } else {
        (second, second_from)
    }
}

fn start(tables: &Tables, mode: Mode) -> (usize, usize, u8, i64) {
    let rows = tables.m.len() / tables.width;
    let width = tables.width;
    match mode {
        Mode::Global => {
            let cell = tables.at(rows - 1, width - 1);
            let (score, state) = best3(tables.m[cell], tables.ix[cell], tables.iy[cell]);
            (rows - 1, width - 1, state, score)
        }
        Mode::Local => {
            let mut best = (0usize, 0usize, STOP, 0_i64);
            for i in 0..rows {
                for j in 0..width {
                    let value = tables.m[tables.at(i, j)];
                    if value > best.3 {
                        best = (i, j, FROM_M, value);
                    }
                }
            }
            best
        }
        Mode::SemiGlobal => {
            let mut best = (0usize, 0usize, STOP, NEG);
            for i in 0..rows {
                consider_end(tables, i, width - 1, &mut best);
            }
            for j in 0..width {
                consider_end(tables, rows - 1, j, &mut best);
            }
            best
        }
    }
}

fn consider_end(tables: &Tables, i: usize, j: usize, best: &mut (usize, usize, u8, i64)) {
    let cell = tables.at(i, j);
    let (value, state) = best3(tables.m[cell], tables.ix[cell], tables.iy[cell]);
    if value > best.3 {
        *best = (i, j, state, value);
    }
}

fn trace(tables: &Tables, mode: Mode, mut i: usize, mut j: usize, mut state: u8) -> Vec<Column> {
    let mut columns = Vec::new();
    while state != STOP && (i > 0 || j > 0) {
        match state {
            FROM_M => {
                if i == 0 || j == 0 {
                    break;
                }
                let from = tables.from_m[tables.at(i, j)];
                columns.push(Column {
                    left: Some(i - 1),
                    right: Some(j - 1),
                });
                i -= 1;
                j -= 1;
                state = from;
            }
            FROM_IX => {
                if i == 0 {
                    break;
                }
                let from = tables.from_ix[tables.at(i, j)];
                columns.push(Column {
                    left: Some(i - 1),
                    right: None,
                });
                i -= 1;
                state = from;
            }
            FROM_IY => {
                if j == 0 {
                    break;
                }
                let from = tables.from_iy[tables.at(i, j)];
                columns.push(Column {
                    left: None,
                    right: Some(j - 1),
                });
                j -= 1;
                state = from;
            }
            _ => break,
        }
        if mode == Mode::SemiGlobal && (i == 0 || j == 0) {
            break;
        }
    }
    columns.reverse();
    columns
}
