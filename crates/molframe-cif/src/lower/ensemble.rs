//! Classifying coordinate models as dense or ragged.
//!
//! Dense storage is used only when every non-coordinate atom-site value agrees
//! row for row. Any difference in identity or annotation requires independent
//! structures so no model's data is silently replaced by the first model's.
//! The `id` serial is the exception: mmCIF files continue it across models
//! (NMR ensembles number atoms 1..N*M), so it is a row counter rather than
//! topology and is not compared. Occupancy and B-factor are per-model values rather than
//! identity, so a difference in them alone stays dense and is reported by
//! [`per_model_values_differ`]; the streaming reader makes the same choice.

use crate::document::Category;

/// One contiguous model in deposition order.
#[derive(Clone, Copy)]
struct ModelRange {
    number: i64,
    start: usize,
    end: usize,
}

/// Deposited model numbers when the category is genuinely ragged.
pub(super) fn ragged_model_numbers(category: &Category, only_first: bool) -> Vec<i64> {
    if only_first {
        return Vec::new();
    }
    let ranges = ranges(category);
    let Some(first) = ranges.first() else {
        return Vec::new();
    };
    if ranges
        .iter()
        .skip(1)
        .all(|candidate| same_atoms(category, *first, *candidate))
    {
        return Vec::new();
    }
    ranges.iter().map(|range| range.number).collect()
}

fn ranges(category: &Category) -> Vec<ModelRange> {
    let mut ranges: Vec<ModelRange> = Vec::new();
    for row in 0..category.row_count() {
        let candidate = category
            .value("pdbx_PDB_model_num", row)
            .and_then(crate::document::CifValue::as_integer);
        let Some(number) = candidate else {
            push_range(&mut ranges, 1, row);
            continue;
        };
        push_range(&mut ranges, number, row);
    }
    ranges
}

fn push_range(ranges: &mut Vec<ModelRange>, number: i64, row: usize) {
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

fn same_atoms(category: &Category, left: ModelRange, right: ModelRange) -> bool {
    if left.end - left.start != right.end - right.start {
        return false;
    }
    let rows = (left.start..left.end).zip(right.start..right.end);
    for (left_row, right_row) in rows {
        for item in category.items() {
            if is_frame_coordinate(item) || item == "id" || is_per_model_value(item) {
                continue;
            }
            if category.value(item, left_row) != category.value(item, right_row) {
                return false;
            }
        }
    }
    true
}

/// Whether the models of a dense category carry different occupancies or B-factors, so the
/// reader must say that model 1's values were kept.
pub(super) fn per_model_values_differ(category: &Category, only_first: bool) -> bool {
    if only_first {
        return false;
    }
    let ranges = ranges(category);
    let Some(first) = ranges.first().copied() else {
        return false;
    };
    ranges.iter().skip(1).any(|candidate| {
        let mut rows = (first.start..first.end).zip(candidate.start..candidate.end);
        rows.any(|(left, right)| {
            category
                .items()
                .filter(|item| is_per_model_value(item))
                .any(|item| category.value(item, left) != category.value(item, right))
        })
    })
}

fn is_per_model_value(item: &str) -> bool {
    matches!(item, "occupancy" | "B_iso_or_equiv")
}

fn is_frame_coordinate(item: &str) -> bool {
    matches!(
        item,
        "Cartn_x" | "Cartn_y" | "Cartn_z" | "pdbx_PDB_model_num"
    )
}

#[cfg(test)]
#[path = "ensemble_tests.rs"]
mod tests;
