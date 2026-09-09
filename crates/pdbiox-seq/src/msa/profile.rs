//! Affine profile-profile alignment and sum-of-pairs scoring.

use super::MsaError;
use crate::{Score, Scoring};

mod memory;
mod workspace;

pub(in crate::msa) use memory::alignment_peak_bytes;
pub(super) use workspace::align;

pub(super) const GAP: u8 = b'-';

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

    pub(super) fn storage_bytes(&self) -> Option<usize> {
        let row_storage = self
            .rows
            .iter()
            .try_fold(0_usize, |total, (_, row)| total.checked_add(row.capacity()))?;
        self.rows
            .capacity()
            .checked_mul(size_of::<(usize, Vec<u8>)>())?
            .checked_add(row_storage)
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
