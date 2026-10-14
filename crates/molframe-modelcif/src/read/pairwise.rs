//! Streams the pair-per-row PAE table straight into a dense matrix.
//!
//! Rows are validated as they arrive and their value is stored in the
//! narrowest exact cell width, so the n-squared compact table is never
//! built. The stream gives up on the first row that is not a lossless
//! fit (several models, gaps, extra items, inexact values, budget) and can
//! then regenerate every row it consumed, so the caller falls back to the
//! compact table without having kept the input.

use super::pairwise_field::Field;
use crate::dense::KNOWN_ITEMS;
use crate::{DensePae, PaeMatrix, PaePolicy, PrecisionProbe, ResidueIndex};
use molframe_cif::CifScalar;
use molframe_core::numeric::f64_to_i64;
use std::collections::HashSet;

const ORDINAL: usize = 0;
const MODEL: usize = 1;
const CHAIN: [usize; 2] = [2, 5];
const SEQUENCE: [usize; 2] = [3, 6];
const COMPONENT: [usize; 2] = [4, 7];
const METRIC: usize = 8;
const VALUE: usize = 9;
const REQUIRED: [usize; 7] = [MODEL, 2, 3, 5, 6, METRIC, VALUE];

type Key<'a> = (&'a str, i64, Option<&'a str>);

/// Incremental dense reader for `ma_qa_metric_local_pairwise`.
#[derive(Debug)]
pub(super) struct PairwiseStream {
    slot: usize,
    policy: PaePolicy,
    order: Vec<usize>,
    closed: bool,
    with_components: bool,
    position: usize,
    pending: [Field; 10],
    stray: Option<(Box<str>, Field)>,
    rejected: bool,
    rows: usize,
    model: Field,
    metric: Field,
    residues: ResidueIndex,
    seen: HashSet<(Box<str>, i64)>,
    first_row: Vec<f64>,
    probe: PrecisionProbe,
    saw_float: bool,
    matrix: Option<PaeMatrix>,
}

impl PairwiseStream {
    /// Starts a stream whose table sits after `slot` other categories.
    pub(super) fn new(slot: usize) -> Self {
        Self {
            slot,
            policy: PaePolicy::new(),
            order: Vec::new(),
            closed: false,
            with_components: false,
            position: 0,
            pending: Default::default(),
            stray: None,
            rejected: false,
            rows: 0,
            model: Field::default(),
            metric: Field::default(),
            residues: ResidueIndex::new(),
            seen: HashSet::new(),
            first_row: Vec::new(),
            probe: PrecisionProbe::new(),
            saw_float: false,
            matrix: None,
        }
    }

    pub(super) const fn slot(&self) -> usize {
        self.slot
    }

    pub(super) fn metric_id(&self) -> Option<String> {
        self.metric.identifier()
    }

    /// Consumes one value; `false` means the table is not a lossless fit.
    pub(super) fn feed(&mut self, item: &str, value: CifScalar<'_>) -> bool {
        let slot = KNOWN_ITEMS.iter().position(|known| *known == item);
        if !self.closed {
            match slot {
                None => return self.reject(item, value),
                Some(slot) if self.order.contains(&slot) => {
                    if self.order.first() != Some(&slot) || !self.finish_first_row() {
                        return self.reject(item, value);
                    }
                }
                Some(slot) => {
                    self.order.push(slot);
                    self.pending[slot].set(value);
                    self.position = self.order.len();
                    return true;
                }
            }
        }
        if slot != self.order.get(self.position).copied() {
            return self.reject(item, value);
        }
        if let Some(slot) = slot {
            self.pending[slot].set(value);
        }
        self.position += 1;
        if self.position == self.order.len() {
            if !self.process_row() {
                self.rejected = true;
                return false;
            }
            self.position = 0;
        }
        true
    }

    /// Whether the whole table was a lossless dense matrix.
    pub(super) fn complete(&mut self) -> bool {
        if self.rejected {
            return false;
        }
        let fits = self.finish_table();
        self.rejected = !fits;
        fits
    }

    /// Re-emits every consumed row, in input order, as scalars.
    pub(super) fn replay<E>(
        &self,
        emit: &mut impl FnMut(&str, CifScalar<'_>) -> Result<(), E>,
    ) -> Result<(), E> {
        for row in 0..self.rows {
            for slot in &self.order {
                emit(KNOWN_ITEMS[*slot], self.cell(row, *slot))?;
            }
        }
        for slot in self.order.iter().take(self.position) {
            emit(KNOWN_ITEMS[*slot], self.pending[*slot].scalar())?;
        }
        if let Some((name, field)) = &self.stray {
            emit(name, field.scalar())?;
        }
        Ok(())
    }

    /// The dense table, once [`complete`](Self::complete) returned `true`.
    pub(super) fn into_dense(self) -> Option<DensePae> {
        let matrix = self.matrix?;
        let items = self
            .order
            .iter()
            .map(|slot| Box::<str>::from(KNOWN_ITEMS[*slot]))
            .collect();
        Some(DensePae {
            model_id: self.model.identifier()?.into(),
            metric_id: self.metric.identifier()?.into(),
            residues: self.residues,
            items,
            slot: self.slot,
            matrix,
        })
    }

    fn reject(&mut self, item: &str, value: CifScalar<'_>) -> bool {
        let mut field = Field::default();
        field.set(value);
        self.stray = Some((item.into(), field));
        self.rejected = true;
        false
    }

    fn finish_first_row(&mut self) -> bool {
        let known = |slot: &usize| self.order.contains(slot);
        let components = [COMPONENT[0], COMPONENT[1]];
        let both = components.iter().all(known);
        if !REQUIRED.iter().all(known) || (!both && components.iter().any(known)) {
            return false;
        }
        self.with_components = both;
        self.closed = true;
        let accepted = self.process_row();
        if accepted {
            self.position = 0;
        }
        accepted
    }

    fn finish_table(&mut self) -> bool {
        if !self.closed {
            if self.order.is_empty() || !self.finish_first_row() {
                return false;
            }
        } else if self.position != 0 {
            return false;
        }
        if self.matrix.is_none() && (self.residues.len() != 1 || !self.start_matrix()) {
            return false;
        }
        self.matrix
            .as_ref()
            .is_some_and(|matrix| matrix.size().checked_mul(matrix.size()) == Some(self.rows))
    }

    fn process_row(&mut self) -> bool {
        let row = self.rows;
        if self.order.contains(&ORDINAL)
            && self.pending[ORDINAL].integer() != i64::try_from(row + 1).ok()
        {
            return false;
        }
        let Some(value) = self.pending[VALUE].number() else {
            return false;
        };
        if row == 0 {
            self.model.copy_from(&self.pending[MODEL]);
            self.metric.copy_from(&self.pending[METRIC]);
        } else if !self.model.same(&self.pending[MODEL]) || !self.metric.same(&self.pending[METRIC])
        {
            return false;
        }
        if self.matrix.is_none() {
            match self.collect_first_row(row, value) {
                Some(true) => {
                    self.saw_float |= self.pending[VALUE].is_float();
                    self.rows += 1;
                    return true;
                }
                Some(false) => {}
                None => return false,
            }
            if !self.start_matrix() {
                return false;
            }
        }
        if !self.place(row, value) {
            return false;
        }
        self.saw_float |= self.pending[VALUE].is_float();
        self.rows += 1;
        true
    }

    /// Row zero of the matrix: `Some(true)` consumed, `Some(false)` row zero ended.
    fn collect_first_row(&mut self, row: usize, value: f64) -> Option<bool> {
        let first = key(&self.pending, 0, self.with_components)?;
        let second = key(&self.pending, 1, self.with_components)?;
        if row == 0 {
            if first != second {
                return None;
            }
        } else if !self.residues.get(0).is_some_and(|known| is(first, &known)) {
            return Some(false);
        }
        if !self.seen.insert((second.0.into(), second.1)) {
            return None;
        }
        let side = self.residues.len().checked_add(1)?;
        if u64::try_from(side.checked_mul(side)?).ok()? > self.policy.max_bytes {
            return None;
        }
        self.residues.push(second.0, second.1, second.2);
        self.first_row.push(value);
        Some(true)
    }

    fn start_matrix(&mut self) -> bool {
        let size = self.residues.len();
        for value in &self.first_row {
            self.probe.observe_f64(*value);
        }
        let Ok(layout) = self.policy.layout(&self.probe, size) else {
            return false;
        };
        let Ok(mut matrix) = PaeMatrix::zeroed(size, layout, self.policy.max_bytes) else {
            return false;
        };
        for (column, value) in self.first_row.iter().enumerate() {
            if matrix.set_f64(0, column, *value).is_err() {
                return false;
            }
        }
        self.first_row = Vec::new();
        self.seen = HashSet::new();
        self.matrix = Some(matrix);
        true
    }

    fn place(&mut self, row: usize, value: f64) -> bool {
        let Some(size) = self.matrix.as_ref().map(PaeMatrix::size) else {
            return false;
        };
        let (line, column) = (row / size, row % size);
        let first = key(&self.pending, 0, self.with_components);
        let second = key(&self.pending, 1, self.with_components);
        let matches =
            |side: Option<Key<'_>>, position: usize| match (side, self.residues.get(position)) {
                (Some(key), Some(known)) => is(key, &known),
                _ => false,
            };
        if line >= size || !matches(first, line) || !matches(second, column) {
            return false;
        }
        self.store(line, column, value)
    }

    fn store(&mut self, line: usize, column: usize, value: f64) -> bool {
        self.probe.observe_f64(value);
        let Some(matrix) = self.matrix.as_mut() else {
            return false;
        };
        if matrix.set_f64(line, column, value).is_err() {
            return false;
        }
        matrix.is_exact() || self.widen(line, column, value)
    }

    /// Re-encodes the cells stored so far in a wider exact layout.
    fn widen(&mut self, line: usize, column: usize, value: f64) -> bool {
        let Some(old) = self.matrix.as_ref() else {
            return false;
        };
        let size = old.size();
        let Ok(layout) = self.policy.layout(&self.probe, size) else {
            return false;
        };
        if layout == old.layout() {
            return false;
        }
        let Ok(mut wider) = PaeMatrix::zeroed(size, layout, self.policy.max_bytes) else {
            return false;
        };
        for index in 0..line * size + column {
            let (row, col) = (index / size, index % size);
            let Some(stored) = old.get_f64(row, col) else {
                return false;
            };
            if wider.set_f64(row, col, stored).is_err() {
                return false;
            }
        }
        if wider.set_f64(line, column, value).is_err() || !wider.is_exact() {
            return false;
        }
        self.matrix = Some(wider);
        true
    }

    fn cell(&self, row: usize, slot: usize) -> CifScalar<'_> {
        let (line, column) = match &self.matrix {
            Some(matrix) => (row / matrix.size(), row % matrix.size()),
            None => (0, row),
        };
        match slot {
            ORDINAL => i64::try_from(row + 1).map_or(CifScalar::Unknown, CifScalar::Integer),
            MODEL => self.model.scalar(),
            METRIC => self.metric.scalar(),
            VALUE => self.value(line, column, row),
            _ => {
                let Some(known) = self.residues.get(if slot > 4 { column } else { line }) else {
                    return CifScalar::Unknown;
                };
                if CHAIN.contains(&slot) {
                    CifScalar::Text(known.chain)
                } else if SEQUENCE.contains(&slot) {
                    CifScalar::Integer(known.sequence)
                } else {
                    known.component.map_or(CifScalar::Unknown, CifScalar::Text)
                }
            }
        }
    }

    fn value(&self, line: usize, column: usize, row: usize) -> CifScalar<'_> {
        let stored = match &self.matrix {
            Some(matrix) => matrix.get_f64(line, column),
            None => self.first_row.get(row).copied(),
        };
        match stored {
            Some(value) if self.saw_float => CifScalar::Float(value),
            Some(value) => CifScalar::Integer(f64_to_i64(value)),
            None => CifScalar::Unknown,
        }
    }
}

fn key(pending: &[Field; 10], side: usize, components: bool) -> Option<Key<'_>> {
    let component = if components {
        Some(pending[COMPONENT[side]].text()?)
    } else {
        None
    };
    Some((
        pending[CHAIN[side]].text()?,
        pending[SEQUENCE[side]].integer()?,
        component,
    ))
}

fn is(key: Key<'_>, known: &crate::ResidueRef<'_>) -> bool {
    key == (known.chain, known.sequence, known.component)
}

#[cfg(test)]
#[path = "pairwise_tests.rs"]
mod tests;
