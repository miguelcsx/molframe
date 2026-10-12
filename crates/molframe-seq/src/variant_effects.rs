//! Variant scores summed up per position, and the check of a table's numbering against a
//! structure's residues.
//!
//! Positions in a variant table are in the assay's own numbering (a `UniProt` entry, a construct)
//! and a structure's residue numbers are the depositor's: they often differ by an offset. So the
//! mapping is checked, never assumed: at every mapped position the table's wild-type letter must
//! be the structure's residue, and the share that agree is reported with the offset used.

use crate::variants::VariantTable;
use std::collections::BTreeMap;

/// What the single substitutions at one position say, for one column.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PositionEffect {
    /// The position, in the table's numbering.
    pub position: i32,
    /// The wild-type residue the table gives.
    pub wild_type: char,
    /// Single substitutions with a value.
    pub count: usize,
    /// Their mean.
    pub mean: f64,
    /// Their median (the mean of the middle two when the count is even).
    pub median: f64,
    /// The share of them on the deleterious side of the [`Cutoff`] used, 0 to 1.
    pub fraction_deleterious: f64,
    /// Lowest value.
    pub min: f64,
    /// Highest value.
    pub max: f64,
}

/// Where a value turns deleterious. A table's scores have no universal direction (a deep
/// scan's fitness falls with damage, a stability ΔΔG rises), so the threshold and the side are
/// stated, never implied.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Cutoff {
    /// The value at which a substitution counts as deleterious (inclusive).
    pub threshold: f64,
    /// Values at or below the threshold are deleterious (else at or above).
    pub low_is_deleterious: bool,
}

impl Default for Cutoff {
    /// Zero, low side: a scan's scores relative to the wild type (`ProteinGym`'s convention).
    fn default() -> Self {
        Self {
            threshold: 0.0,
            low_is_deleterious: true,
        }
    }
}

impl Cutoff {
    /// Whether `value` is deleterious.
    #[must_use]
    pub fn is_deleterious(&self, value: f64) -> bool {
        if self.low_is_deleterious {
            value <= self.threshold
        } else {
            value >= self.threshold
        }
    }

    /// One sentence naming the rule.
    #[must_use]
    pub fn text(&self) -> String {
        format!(
            "deleterious = value {} {}",
            if self.low_is_deleterious { "<=" } else { ">=" },
            self.threshold
        )
    }
}

/// A count as a float; counts past `u32::MAX` saturate, far beyond any table.
fn count_as_f64(count: usize) -> f64 {
    f64::from(u32::try_from(count).unwrap_or(u32::MAX))
}

/// The median of an ascending slice (the mean of the middle two when the length is even), NaN
/// when it is empty.
#[must_use]
pub fn median_of_sorted(sorted: &[f64]) -> f64 {
    let n = sorted.len();
    match (n, sorted.get(n / 2)) {
        (0, _) | (_, None) => f64::NAN,
        (_, Some(upper)) if n % 2 == 1 => *upper,
        (_, Some(upper)) => sorted
            .get(n / 2 - 1)
            .map_or(f64::NAN, |lower| (lower + upper) / 2.0),
    }
}

/// [`per_position_at`] with the default [`Cutoff`].
#[must_use]
pub fn per_position(table: &VariantTable, column: usize) -> Vec<PositionEffect> {
    per_position_at(table, column, Cutoff::default())
}

/// The single substitutions of column `column`, summed up by position (variants of several
/// substitutions are not attributable to one position and are left out of the sums).
#[must_use]
pub fn per_position_at(table: &VariantTable, column: usize, cutoff: Cutoff) -> Vec<PositionEffect> {
    let mut by: BTreeMap<i32, (char, Vec<f64>)> = BTreeMap::new();
    for variant in &table.variants {
        let ([change], Some(Some(value))) = (variant.changes.as_slice(), variant.values.get(column))
        else {
            continue;
        };
        by.entry(change.position)
            .or_insert((change.wild_type, Vec::new()))
            .1
            .push(*value);
    }
    by.into_iter()
        .filter_map(|(position, (wild_type, mut values))| {
            let bad = values.iter().filter(|v| cutoff.is_deleterious(**v)).count();
            values.sort_by(f64::total_cmp);
            let (min, max) = (*values.first()?, *values.last()?);
            let count = count_as_f64(values.len());
            Some(PositionEffect {
                position,
                wild_type,
                count: values.len(),
                mean: values.iter().sum::<f64>() / count,
                median: median_of_sorted(&values),
                fraction_deleterious: count_as_f64(bad) / count,
                min,
                max,
            })
        })
        .collect()
}

/// How the table's positions fell on a chain's residues.
#[derive(Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Mapping {
    /// Added to the table's positions to get the structure's residue numbers.
    pub offset: i32,
    /// Positions whose wild-type letter is the structure's residue.
    pub agree: usize,
    /// Positions whose letter is another residue (the numbering is wrong or the construct differs).
    pub disagree: usize,
    /// Positions the chain has no residue for.
    pub absent: usize,
}

impl Mapping {
    /// The share of mapped positions whose letters agree, 0 to 1.
    #[must_use]
    pub fn agreement(&self) -> Option<f64> {
        let mapped = self.agree + self.disagree;
        (mapped > 0).then(|| count_as_f64(self.agree) / count_as_f64(mapped))
    }
}

/// Checks `effects` against the residues `(number, letter)` of a chain at `offset`.
#[must_use]
pub fn map_positions(effects: &[PositionEffect], residues: &[(i32, char)], offset: i32) -> Mapping {
    let by_number: BTreeMap<i32, char> = residues.iter().copied().collect();
    let mut mapping = Mapping {
        offset,
        agree: 0,
        disagree: 0,
        absent: 0,
    };
    for effect in effects {
        match by_number.get(&effect.position.saturating_add(offset)) {
            None => mapping.absent += 1,
            Some(c) if *c == effect.wild_type => mapping.agree += 1,
            Some(_) => mapping.disagree += 1,
        }
    }
    mapping
}

/// The offset in `-range..=range` (residues either way) at which the wild-type letters agree
/// best, with its mapping. Ties go to the smallest offset in magnitude; `None` when nothing
/// agrees anywhere.
#[must_use]
pub fn best_offset(
    effects: &[PositionEffect],
    residues: &[(i32, char)],
    range: i32,
) -> Option<Mapping> {
    let mut best: Option<Mapping> = None;
    let range = range.max(0);
    // Nearest first (0, -1, 1, -2, ...), so a tie keeps the smaller magnitude.
    let offsets = (0..=range).flat_map(|k| [-k, k].into_iter().take(if k == 0 { 1 } else { 2 }));
    for offset in offsets {
        let mapping = map_positions(effects, residues, offset);
        if mapping.agree == 0 {
            continue;
        }
        let better = best.as_ref().is_none_or(|b| {
            mapping.agree > b.agree
                || (mapping.agree == b.agree && mapping.offset.abs() < b.offset.abs())
        });
        if better {
            best = Some(mapping);
        }
    }
    best
}

#[cfg(test)]
#[path = "variant_effects_tests.rs"]
mod tests;
