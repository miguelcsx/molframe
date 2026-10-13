//! Conversion of a compact pair-per-row PAE table into a dense matrix.
//!
//! The table is converted only when nothing is lost: one model, one PAE
//! metric, every ordered pair stated once in row-major order, only the
//! dictionary's pairwise items, consecutive ordinals, and values the chosen
//! cell width reproduces exactly. Anything else stays compact and is written
//! back unchanged.

use super::error::ConfidenceError;
use super::layout::PaePolicy;
use super::pae::PaeMatrix;
use super::pairs::DensePae;
use super::precision::PrecisionProbe;
use super::residues::ResidueIndex;
use crate::{CompactColumn, ModelCategory, ModelCif};

pub(crate) const PAIRWISE: &str = "ma_qa_metric_local_pairwise";

pub(crate) const KNOWN_ITEMS: [&str; 10] = [
    "ordinal_id",
    "model_id",
    "label_asym_id_1",
    "label_seq_id_1",
    "label_comp_id_1",
    "label_asym_id_2",
    "label_seq_id_2",
    "label_comp_id_2",
    "metric_id",
    "metric_value",
];

struct Columns<'a> {
    ordinal: Option<&'a CompactColumn>,
    model: &'a CompactColumn,
    metric: &'a CompactColumn,
    value: &'a CompactColumn,
    chains: [&'a CompactColumn; 2],
    sequences: [&'a CompactColumn; 2],
    components: [Option<&'a CompactColumn>; 2],
}

type Key<'a> = (&'a str, i64, Option<&'a str>);

impl<'a> Columns<'a> {
    fn new(category: &'a ModelCategory) -> Option<Self> {
        if category
            .items()
            .iter()
            .any(|item| !KNOWN_ITEMS.contains(&&**item))
        {
            return None;
        }
        Some(Self {
            ordinal: category.column("ordinal_id"),
            model: category.column("model_id")?,
            metric: category.column("metric_id")?,
            value: category.column("metric_value")?,
            chains: [
                category.column("label_asym_id_1")?,
                category.column("label_asym_id_2")?,
            ],
            sequences: [
                category.column("label_seq_id_1")?,
                category.column("label_seq_id_2")?,
            ],
            components: [
                category.column("label_comp_id_1"),
                category.column("label_comp_id_2"),
            ],
        })
    }

    fn residue(&self, side: usize, row: usize) -> Option<Key<'a>> {
        let chain = self.chains[side].identifier(row)?;
        let sequence = self.sequences[side].value(row)?.integer()?;
        let component = match self.components[side] {
            Some(column) => Some(column.identifier(row)?),
            None => None,
        };
        Some((chain, sequence, component))
    }

    fn number(&self, row: usize) -> Option<f64> {
        self.value
            .value(row)?
            .number()
            .filter(|value| value.is_finite())
    }

    fn ordinal_in_order(&self, row: usize) -> bool {
        self.ordinal.is_none_or(|column| {
            column
                .value(row)
                .and_then(crate::ValueRef::integer)
                .and_then(|value| usize::try_from(value).ok())
                == Some(row + 1)
        })
    }
}

/// What the first pass learns about a convertible table.
struct Survey<'a> {
    model_id: &'a str,
    metric_id: &'a str,
    residues: ResidueIndex,
    probe: PrecisionProbe,
}

/// Integer square root when `rows` is a perfect square.
fn side(rows: usize) -> Option<usize> {
    let mut size = 0_usize;
    while size.checked_mul(size)? < rows {
        size += 1;
    }
    (size * size == rows).then_some(size)
}

/// First pass: one model and metric, residue order, value precision, row order.
fn survey<'a>(columns: &Columns<'a>, rows: usize, model: &ModelCif) -> Option<Survey<'a>> {
    let size = side(rows).filter(|size| *size > 0)?;
    let model_id = columns.model.identifier(0)?;
    let metric_id = columns.metric.identifier(0)?;
    let is_pae = model.confidence().definitions().any(|definition| {
        definition.id == metric_id && definition.metric_type.eq_ignore_ascii_case("PAE")
    });
    if !is_pae {
        return None;
    }
    let mut residues = ResidueIndex::new();
    for row in 0..size {
        let (chain, sequence, component) = columns.residue(1, row)?;
        if residues.position(chain, sequence).is_some() {
            return None;
        }
        residues.push(chain, sequence, component);
    }
    let mut probe = PrecisionProbe::new();
    for row in 0..rows {
        if columns.model.identifier(row)? != model_id
            || columns.metric.identifier(row)? != metric_id
            || !columns.ordinal_in_order(row)
            || !matches(&residues, row / size, columns.residue(0, row)?)
            || !matches(&residues, row % size, columns.residue(1, row)?)
        {
            return None;
        }
        probe.observe_f64(columns.number(row)?);
    }
    Some(Survey {
        model_id,
        metric_id,
        residues,
        probe,
    })
}

fn matches(residues: &ResidueIndex, position: usize, key: Key<'_>) -> bool {
    residues
        .get(position)
        .is_some_and(|known| (known.chain, known.sequence, known.component) == key)
}

impl ModelCif {
    /// Replaces a convertible pair-per-row PAE table with a dense matrix.
    ///
    /// Returns whether the table was converted. A table that cannot be
    /// represented exactly (several models, gaps, extra items, values the cell
    /// width cannot hold) is left compact and is not an error.
    ///
    /// # Errors
    ///
    /// Returns the budget or allocation failure of an otherwise convertible table.
    pub fn densify_pae(&mut self, policy: PaePolicy) -> Result<bool, ConfidenceError> {
        let Some((slot, dense)) = build(self, policy)? else {
            return Ok(false);
        };
        self.categories.remove(slot);
        self.pae = Some(dense);
        Ok(true)
    }

    /// Swaps the stored matrix for one of the same side, such as a requantised copy.
    ///
    /// The matrix is moved in, not copied. Writers then report its quantisation
    /// error when it is not exact.
    ///
    /// # Errors
    ///
    /// Returns a size mismatch when there is no dense PAE or the sides differ.
    pub fn replace_pae_matrix(&mut self, matrix: PaeMatrix) -> Result<(), ConfidenceError> {
        let current = self.pae.as_ref().map_or(0, |dense| dense.matrix.size());
        match self.pae.as_mut() {
            Some(dense) if current == matrix.size() => {
                dense.matrix = matrix;
                Ok(())
            }
            _ => Err(ConfidenceError::SizeMismatch {
                residues: current,
                matrix: matrix.size(),
            }),
        }
    }
}

fn build(
    model: &ModelCif,
    policy: PaePolicy,
) -> Result<Option<(usize, DensePae)>, ConfidenceError> {
    let Some(slot) = model
        .categories
        .iter()
        .position(|category| category.name() == PAIRWISE)
    else {
        return Ok(None);
    };
    let category = &model.categories[slot];
    let Some(columns) = Columns::new(category) else {
        return Ok(None);
    };
    let Some(shape) = survey(&columns, category.row_count(), model) else {
        return Ok(None);
    };
    let size = shape.residues.len();
    let layout = policy.layout(&shape.probe, size)?;
    let mut matrix = PaeMatrix::zeroed(size, layout, policy.max_bytes)?;
    for row in 0..category.row_count() {
        let Some(value) = columns.number(row) else {
            return Ok(None);
        };
        matrix.set_f64(row / size, row % size, value)?;
    }
    if !matrix.is_exact() {
        return Ok(None);
    }
    let dense = DensePae {
        model_id: shape.model_id.into(),
        metric_id: shape.metric_id.into(),
        residues: shape.residues,
        items: category.items().to_vec(),
        slot,
        matrix,
    };
    Ok(Some((slot, dense)))
}

#[cfg(test)]
#[path = "densify_tests.rs"]
mod tests;
