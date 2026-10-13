//! Shared Gotoh recurrence with rolling score rows and a packed traceback.
//!
//! Scores live in two rolling rows. The only table kept for the whole
//! alignment is one byte per cell holding the three predecessor pointers, so
//! a local or semi-global alignment costs one byte per cell instead of
//! twenty-seven.

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

/// Three 2-bit predecessor pointers per cell, `M | Ix << 2 | Iy << 4`.
const ALL_STOP: u8 = STOP | (STOP << 2) | (STOP << 4);

struct Trace {
    rows: Vec<Row>,
    bits: Vec<u8>,
}

struct Row {
    offset: usize,
    first: usize,
    end: usize,
}

impl Trace {
    fn at(&self, row: usize, col: usize) -> usize {
        let layout = &self.rows[row];
        if col == 0 {
            return layout.offset;
        }
        if col < layout.first || col >= layout.end {
            return 0;
        }
        layout.offset + 1 + col - layout.first
    }

    fn pointer(&self, row: usize, col: usize, shift: u8) -> u8 {
        (self.bits[self.at(row, col)] >> shift) & 3
    }

    fn set(&mut self, row: usize, col: usize, from: [u8; 3]) {
        let index = self.at(row, col);
        self.bits[index] = from[0] | (from[1] << 2) | (from[2] << 4);
    }
}

/// One DP row of the three Gotoh states; `lo..hi` is the span written last.
struct Line {
    m: Vec<i64>,
    ix: Vec<i64>,
    iy: Vec<i64>,
    lo: usize,
    hi: usize,
}

impl Line {
    fn new(width: usize) -> Self {
        Self {
            m: vec![NEG; width],
            ix: vec![NEG; width],
            iy: vec![NEG; width],
            lo: 0,
            hi: 0,
        }
    }

    fn clear(&mut self) {
        for column in self.lo..self.hi {
            self.m[column] = NEG;
            self.ix[column] = NEG;
            self.iy[column] = NEG;
        }
    }

    fn best(&self, column: usize) -> (i64, u8) {
        best3(self.m[column], self.ix[column], self.iy[column])
    }
}

type End = (usize, usize, u8, i64);

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
    let mut trace = layout(rows, width, band)?;
    let cell = CellParameters {
        left,
        right,
        scorer,
        gap_open,
        gap_extend,
        mode,
    };
    let mut previous = Line::new(width);
    let mut current = Line::new(width);
    seed_row_zero(&mut previous, width, &cell)?;
    seed_boundary_pointers(&mut trace, rows, width, mode);

    let mut best_local: End = (0, 0, STOP, 0);
    let mut best_column: End = (0, 0, STOP, NEG);
    if mode == Mode::SemiGlobal {
        consider(&previous, 0, width - 1, &mut best_column);
    }

    // A band restricts which cells are filled, so it must restrict which cells
    // are visited, not skip them inside a full traversal.
    for i in 1..rows {
        let (first, last) = banded_columns(i, width, band);
        current.clear();
        seed_column_zero(&mut current, i, &cell)?;
        current.lo = 0;
        current.hi = last.clamp(1, width);
        for j in first..last {
            fill_cell(&mut trace, &previous, &mut current, &cell, i, j)?;
            if mode == Mode::Local && current.m[j] > best_local.3 {
                best_local = (i, j, FROM_M, current.m[j]);
            }
        }
        if mode == Mode::SemiGlobal {
            consider(&current, i, width - 1, &mut best_column);
        }
        std::mem::swap(&mut previous, &mut current);
    }

    let (start_row, start_col, start_state, score) = match mode {
        Mode::Global => {
            let (score, state) = previous.best(width - 1);
            (rows - 1, width - 1, state, score)
        }
        Mode::Local => best_local,
        Mode::SemiGlobal => {
            for j in 0..width {
                consider(&previous, rows - 1, j, &mut best_column);
            }
            best_column
        }
    };
    if score <= NEG {
        return Err(AlignError::NoAlignmentPath);
    }
    let columns = trace_back(&trace, mode, start_row, start_col, start_state);
    Ok(Alignment {
        score: i32::try_from(score).map_err(|_| AlignError::NumericOverflow)?,
        columns,
    })
}

fn consider(line: &Line, row: usize, column: usize, best: &mut End) {
    let (value, state) = line.best(column);
    if value > best.3 {
        *best = (row, column, state, value);
    }
}

/// The half-open column range row `i` fills under an optional band.
///
/// Without a band that is every column. With one it is the columns within
/// `band` of the diagonal, clamped to the row.
fn banded_columns(row: usize, width: usize, band: Option<usize>) -> (usize, usize) {
    let Some(band) = band else {
        return (1, width);
    };
    let first = row.saturating_sub(band).max(1);
    let last = match row.checked_add(band).and_then(|last| last.checked_add(1)) {
        Some(last) => last.min(width),
        None => width,
    };
    (first, last.max(first))
}

fn layout(rows: usize, width: usize, band: Option<usize>) -> Result<Trace, AlignError> {
    // Slot zero is the immutable unreachable predecessor outside the band.
    // Keep the complete boundary row/column to preserve free-end semantics.
    let mut spans = Vec::with_capacity(rows);
    let mut size = 1usize;
    for row in 0..rows {
        let (first, end) = if row == 0 {
            (1, width)
        } else {
            banded_columns(row, width, band)
        };
        spans.push(Row {
            offset: size,
            first,
            end,
        });
        size = size
            .checked_add(1)
            .and_then(|value| value.checked_add(end - first))
            .ok_or(AlignError::NumericOverflow)?;
    }
    let mut bits = Vec::new();
    bits.try_reserve_exact(size)
        .map_err(|_| AlignError::NumericOverflow)?;
    bits.resize(size, ALL_STOP);
    Ok(Trace { rows: spans, bits })
}

/// Boundary gaps of a global alignment trace back along the edge to the origin.
fn seed_boundary_pointers(trace: &mut Trace, rows: usize, width: usize, mode: Mode) {
    if mode != Mode::Global {
        return;
    }
    for i in 1..rows {
        let from = if i == 1 { FROM_M } else { FROM_IX };
        trace.set(i, 0, [STOP, from, STOP]);
    }
    for j in 1..width {
        let from = if j == 1 { FROM_M } else { FROM_IY };
        trace.set(0, j, [STOP, STOP, from]);
    }
}

fn seed_row_zero<S: Score>(
    line: &mut Line,
    width: usize,
    cell: &CellParameters<'_, S>,
) -> Result<(), AlignError> {
    line.lo = 0;
    line.hi = width;
    line.m[0] = 0;
    match cell.mode {
        Mode::Global => {
            let mut penalty = i64::from(cell.gap_open);
            for j in 1..width {
                line.iy[j] = penalty;
                if j + 1 < width {
                    penalty = penalty
                        .checked_add(i64::from(cell.gap_extend))
                        .ok_or(AlignError::NumericOverflow)?;
                }
            }
        }
        Mode::Local | Mode::SemiGlobal => {
            line.ix[0] = 0;
            line.iy[0] = 0;
            for j in 1..width {
                line.m[j] = 0;
                line.iy[j] = 0;
            }
        }
    }
    Ok(())
}

fn seed_column_zero<S: Score>(
    line: &mut Line,
    row: usize,
    cell: &CellParameters<'_, S>,
) -> Result<(), AlignError> {
    match cell.mode {
        Mode::Global => {
            let steps = i64::try_from(row - 1).map_err(|_| AlignError::NumericOverflow)?;
            line.ix[0] = i64::from(cell.gap_extend)
                .checked_mul(steps)
                .and_then(|extra| extra.checked_add(i64::from(cell.gap_open)))
                .ok_or(AlignError::NumericOverflow)?;
        }
        Mode::Local | Mode::SemiGlobal => {
            line.m[0] = 0;
            line.ix[0] = 0;
        }
    }
    Ok(())
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
    trace: &mut Trace,
    previous: &Line,
    current: &mut Line,
    parameters: &CellParameters<'_, S>,
    i: usize,
    j: usize,
) -> Result<(), AlignError> {
    let (diagonal_best, diagonal_from) = previous.best(j - 1);
    let substitution = parameters
        .scorer
        .score(parameters.left[i - 1], parameters.right[j - 1]);
    let mut match_value = add(diagonal_best, substitution)?;
    let mut match_from = diagonal_from;
    if parameters.mode == Mode::Local && match_value < 0 {
        match_value = 0;
        match_from = STOP;
    }
    current.m[j] = match_value;

    let (vertical_gap_value, vertical_gap_from) = better(
        add(previous.m[j], parameters.gap_open)?,
        FROM_M,
        add(previous.ix[j], parameters.gap_extend)?,
        FROM_IX,
    );
    current.ix[j] = vertical_gap_value;

    let (horizontal_gap_value, horizontal_gap_from) = better(
        add(current.m[j - 1], parameters.gap_open)?,
        FROM_M,
        add(current.iy[j - 1], parameters.gap_extend)?,
        FROM_IY,
    );
    current.iy[j] = horizontal_gap_value;
    trace.set(i, j, [match_from, vertical_gap_from, horizontal_gap_from]);
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

fn trace_back(
    tables: &Trace,
    mode: Mode,
    mut i: usize,
    mut j: usize,
    mut state: u8,
) -> Vec<Column> {
    let mut columns = Vec::new();
    while state != STOP && (i > 0 || j > 0) {
        match state {
            FROM_M => {
                if i == 0 || j == 0 {
                    break;
                }
                let from = tables.pointer(i, j, 0);
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
                let from = tables.pointer(i, j, 2);
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
                let from = tables.pointer(i, j, 4);
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

#[cfg(test)]
#[path = "dynamic_tests.rs"]
mod tests;
