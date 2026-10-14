//! Dense-versus-ragged classification over columnar values.
//!
//! A multi-model file is dense only when every non-coordinate value matches
//! the first model row for row. Occupancy and B-factor are per-model values, not
//! topology: they may differ without making the file ragged, and the layout says
//! so for the caller to report. The `id` serial is numbering that mmCIF continues
//! across models, so it is not compared at all. Extra `atom_site` columns are decoded one at a
//! time and dropped immediately, so classification does not create a DOM.

use super::column::{AtomColumns, DecodedColumn, is_frame_coordinate};
use super::container::DirectColumn;
use molframe_core::diagnostic::Diagnostic;

#[derive(Clone, Copy)]
struct ModelRange {
    number: i64,
    start: usize,
    end: usize,
}

pub(super) struct EnsembleLayout {
    pub(super) model_count: usize,
    pub(super) ragged: bool,
    /// Occupancies or B-factors differ between models of a dense layout.
    pub(super) values_differ: bool,
}

pub(super) fn classify(
    atoms: &AtomColumns,
    extra_columns: Vec<DirectColumn<'_>>,
) -> Result<EnsembleLayout, Diagnostic> {
    let ranges = ranges(atoms);
    if ranges.len() < 2 {
        return Ok(EnsembleLayout {
            model_count: ranges.len(),
            ragged: false,
            values_differ: false,
        });
    }
    if ranges
        .iter()
        .skip(1)
        .any(|range| range.end - range.start != ranges[0].end - ranges[0].start)
    {
        return Ok(EnsembleLayout {
            model_count: ranges.len(),
            ragged: true,
            values_differ: false,
        });
    }
    let mut values_differ = false;
    for (name, column) in atoms.iter() {
        if is_frame_coordinate(name) || name == "id" {
            continue;
        }
        if is_model_value(name) {
            values_differ |= !same_models(column, &ranges);
        } else if !same_models(column, &ranges) {
            return Ok(EnsembleLayout {
                model_count: ranges.len(),
                ragged: true,
                values_differ: false,
            });
        }
    }
    for encoded in extra_columns {
        if is_frame_coordinate(&encoded.name) {
            continue;
        }
        let column = DecodedColumn::decode(&encoded, atoms.row_count(), false)?;
        if !same_models(&column, &ranges) {
            return Ok(EnsembleLayout {
                model_count: ranges.len(),
                ragged: true,
                values_differ: false,
            });
        }
    }
    Ok(EnsembleLayout {
        model_count: ranges.len(),
        ragged: false,
        values_differ,
    })
}

/// Columns that describe a model's values rather than its atoms.
fn is_model_value(name: &str) -> bool {
    matches!(name, "occupancy" | "B_iso_or_equiv")
}

fn ranges(atoms: &AtomColumns) -> Vec<ModelRange> {
    let model = atoms.column("pdbx_PDB_model_num");
    let mut ranges: Vec<ModelRange> = Vec::new();
    for row in 0..atoms.row_count() {
        let number = model_or_one(model.and_then(|column| column.integer(row)));
        if let Some(current) = ranges.last_mut()
            && current.number == number
        {
            current.end = row + 1;
        } else {
            ranges.push(ModelRange {
                number,
                start: row,
                end: row + 1,
            });
        }
    }
    ranges
}

fn model_or_one(number: Option<i64>) -> i64 {
    let Some(number) = number else {
        return 1;
    };
    number
}

fn same_models(column: &DecodedColumn, ranges: &[ModelRange]) -> bool {
    let Some(first) = ranges.first() else {
        return true;
    };
    for candidate in ranges.iter().skip(1) {
        let mut rows = (first.start..first.end).zip(candidate.start..candidate.end);
        if rows.any(|(left, right)| !column.same_rows(left, right)) {
            return false;
        }
    }
    true
}
