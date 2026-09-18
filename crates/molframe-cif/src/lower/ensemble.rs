//! Classifying coordinate models as dense or ragged.
//!
//! Dense storage is used only when every non-coordinate atom-site value agrees
//! row for row. Any difference in identity or annotation requires independent
//! structures so no model's data is silently replaced by the first model's.

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
        let number = match category
            .value("pdbx_PDB_model_num", row)
            .and_then(crate::document::CifValue::as_integer)
        {
            Some(number) => number,
            None => 1,
        };
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

fn same_atoms(category: &Category, left: ModelRange, right: ModelRange) -> bool {
    if left.end - left.start != right.end - right.start {
        return false;
    }
    let rows = (left.start..left.end).zip(right.start..right.end);
    for (left_row, right_row) in rows {
        for item in category.items() {
            if is_frame_coordinate(item) {
                continue;
            }
            if category.value(item, left_row) != category.value(item, right_row) {
                return false;
            }
        }
    }
    true
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
