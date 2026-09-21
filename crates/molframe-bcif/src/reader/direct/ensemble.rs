//! Dense-versus-ragged classification over columnar values.
//!
//! A multi-model file is dense only when every non-coordinate value matches
//! the first model row for row. Extra `atom_site` columns are decoded one at a
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
        });
    }
    for (name, column) in atoms.iter() {
        if !is_frame_coordinate(name) && !same_models(column, &ranges) {
            return Ok(EnsembleLayout {
                model_count: ranges.len(),
                ragged: true,
            });
        }
    }
    for encoded in extra_columns {
        if is_frame_coordinate(&encoded.name) {
            continue;
        }
        let column = DecodedColumn::decode(encoded, atoms.row_count(), false)?;
        if !same_models(&column, &ranges) {
            return Ok(EnsembleLayout {
                model_count: ranges.len(),
                ragged: true,
            });
        }
    }
    Ok(EnsembleLayout {
        model_count: ranges.len(),
        ragged: false,
    })
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
