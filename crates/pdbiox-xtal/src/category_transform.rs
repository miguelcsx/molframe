//! Shared lowering of recorded Cartesian matrix/vector categories.

use crate::AffineTransform;
use pdbiox_cif::Category;
use pdbiox_core::{Code, Diagnostic};

pub(crate) fn read_affine(
    category: &Category,
    row: usize,
    code: Code,
) -> Result<AffineTransform, Diagnostic> {
    let mut matrix = [[0.0; 3]; 3];
    let mut translation = [0.0; 3];
    for (matrix_row, values) in matrix.iter_mut().enumerate() {
        for (column, value) in values.iter_mut().enumerate() {
            let item = format!("matrix[{}][{}]", matrix_row + 1, column + 1);
            *value = required_float(category, &item, row, code)?;
        }
        let item = format!("vector[{}]", matrix_row + 1);
        translation[matrix_row] = required_float(category, &item, row, code)?;
    }
    let transform = AffineTransform::new(matrix, translation);
    if transform.is_finite() {
        Ok(transform)
    } else {
        Err(malformed(category, row, "matrix", code))
    }
}

fn required_float(
    category: &Category,
    item: &str,
    row: usize,
    code: Code,
) -> Result<f64, Diagnostic> {
    category
        .value(item, row)
        .and_then(pdbiox_cif::CifValue::as_float)
        .filter(|value| value.is_finite())
        .ok_or_else(|| malformed(category, row, item, code))
}

pub(crate) fn malformed(category: &Category, row: usize, item: &str, code: Code) -> Diagnostic {
    Diagnostic::new(code)
        .with_context("category", category.name())
        .with_context("item", item)
        .with_context("row", row.to_string())
}
