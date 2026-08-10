//! Linear iteration over encoded columns.

use super::super::bits::unpack_one;
use super::model::{ColumnValue, EncodedColumn};

/// Sequential reader over an encoded column.
#[derive(Debug)]
pub struct ColumnIter<'a, T> {
    column: &'a EncodedColumn<T>,
    position: u32,
    running: i64,
    run: usize,
}

impl<'a, T> ColumnIter<'a, T> {
    pub(super) const fn new(column: &'a EncodedColumn<T>) -> Self {
        Self {
            column,
            position: 0,
            running: 0,
            run: 0,
        }
    }
}

impl<T: ColumnValue> Iterator for ColumnIter<'_, T> {
    type Item = T;

    fn next(&mut self) -> Option<T> {
        if self.position as usize >= self.column.len() {
            return None;
        }
        let value = match self.column {
            EncodedColumn::Plain(values) => copied_at(values, self.position)?,
            EncodedColumn::Constant { value, .. } => *value,
            EncodedColumn::RunLength { values, run_ends } => {
                next_run_length(values, run_ends, self.position, &mut self.run)?
            }
            EncodedColumn::BitPacked { data, width, .. } => {
                T::from_bits(unpack_one(data, *width, self.position)?)?
            }
            EncodedColumn::Delta { first, deltas } => {
                next_delta(*first, deltas, self.position, &mut self.running)?
            }
        };
        self.position = self.position.checked_add(1)?;
        Some(value)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.column.len() - self.position as usize;
        (remaining, Some(remaining))
    }
}

fn copied_at<T: Copy>(values: &[T], position: u32) -> Option<T> {
    values.get(usize::try_from(position).ok()?).copied()
}

fn next_run_length<T: Copy>(
    values: &[T],
    run_ends: &[u32],
    position: u32,
    run: &mut usize,
) -> Option<T> {
    while run_ends.get(*run).is_some_and(|end| *end <= position) {
        *run += 1;
    }
    values.get(*run).copied()
}

fn next_delta<T: ColumnValue>(
    first: T,
    deltas: &[i32],
    position: u32,
    running: &mut i64,
) -> Option<T> {
    if position == 0 {
        *running = first.to_bits().cast_signed();
    } else {
        let index = usize::try_from(position.checked_sub(1)?).ok()?;
        *running += i64::from(*deltas.get(index)?);
    }
    T::from_bits((*running).cast_unsigned())
}

impl<T: ColumnValue> ExactSizeIterator for ColumnIter<'_, T> {}

impl<'a, T: ColumnValue> IntoIterator for &'a EncodedColumn<T> {
    type Item = T;
    type IntoIter = ColumnIter<'a, T>;
    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}
