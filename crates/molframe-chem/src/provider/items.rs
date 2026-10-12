//! Item-level readers shared by the CCD lowering.

use molframe_cif::{Category, CifValue};
use molframe_core::diagnostic::{Code, Diagnostic};

type LoweringResult<T> = Result<T, Vec<Diagnostic>>;

pub(super) fn collect_rows<T>(
    category: &Category,
    parse: impl Fn(usize) -> Result<T, Diagnostic>,
) -> LoweringResult<Vec<T>> {
    let mut values = Vec::with_capacity(category.row_count());
    let mut findings = Vec::new();
    for row in 0..category.row_count() {
        match parse(row) {
            Ok(value) => values.push(value),
            Err(finding) => findings.push(finding),
        }
    }
    if findings.is_empty() {
        Ok(values)
    } else {
        Err(findings)
    }
}

pub(super) fn required_identifier(
    category: &Category,
    item: &str,
    row: usize,
) -> Result<Box<str>, Diagnostic> {
    category
        .identifier(item, row)
        .map(|value| value.into_owned().into())
        .ok_or_else(|| missing_item(category, item, row))
}

pub(super) fn required_text(
    category: &Category,
    item: &str,
    row: usize,
) -> Result<Box<str>, Diagnostic> {
    category
        .text(item, row)
        .map(Into::into)
        .ok_or_else(|| missing_item(category, item, row))
}

pub(super) fn required_integer(
    category: &Category,
    item: &str,
    row: usize,
) -> Result<i64, Diagnostic> {
    let value = category
        .value(item, row)
        .ok_or_else(|| missing_item(category, item, row))?;
    if matches!(value, CifValue::Unknown | CifValue::Inapplicable) {
        return Err(missing_item(category, item, row));
    }
    value
        .as_integer()
        .ok_or_else(|| wrong_type(category, item, row))
}

pub(super) fn required_number(
    category: &Category,
    item: &str,
    row: usize,
) -> Result<f64, Diagnostic> {
    let value = category
        .value(item, row)
        .ok_or_else(|| missing_item(category, item, row))?;
    if matches!(value, CifValue::Unknown | CifValue::Inapplicable) {
        return Err(missing_item(category, item, row));
    }
    value
        .as_float()
        .filter(|number| number.is_finite())
        .ok_or_else(|| wrong_type(category, item, row))
}

pub(super) fn one_finding(finding: Diagnostic) -> Vec<Diagnostic> {
    vec![finding]
}

pub(super) fn required_flag(
    category: &Category,
    item: &str,
    row: usize,
) -> Result<bool, Diagnostic> {
    match category.text(item, row) {
        Some(value) if value.eq_ignore_ascii_case("Y") => Ok(true),
        Some(value) if value.eq_ignore_ascii_case("N") => Ok(false),
        Some(_) => Err(invalid_value(category, item, row)),
        None => Err(missing_item(category, item, row)),
    }
}

pub(super) fn missing_item(category: &Category, item: &str, row: usize) -> Diagnostic {
    Diagnostic::new(Code::E2002)
        .in_category(category.name())
        .with_context("item", item)
        .with_context("row", row.to_string())
}

pub(super) fn wrong_type(category: &Category, item: &str, row: usize) -> Diagnostic {
    Diagnostic::new(Code::E2004)
        .in_category(category.name())
        .with_context("item", item)
        .with_context("row", row.to_string())
}

pub(super) fn invalid_value(category: &Category, item: &str, row: usize) -> Diagnostic {
    Diagnostic::new(Code::E2003)
        .in_category(category.name())
        .with_context("item", item)
        .with_context("row", row.to_string())
}

/// Whether the file wrote `?` (unknown) or `.` (inapplicable) for this item.
///
/// Both are real values in the archive, distinct from a column that is absent.
pub(super) fn unknown_marker(category: &Category, item: &str, row: usize) -> bool {
    matches!(
        category.value(item, row),
        Some(CifValue::Unknown | CifValue::Inapplicable)
    )
}

/// Formal charge of a component atom.
///
/// The dictionary declares `charge` with a default of 0, and the archive writes
/// `?` (unknown) or `.` (inapplicable) for atoms whose charge was never
/// curated, such as the iron of haem. Those two markers therefore resolve to
/// the dictionary default rather than failing the whole component; an
/// explicitly written value is always honoured, and a column that is absent
/// altogether is still a missing item.
pub(super) fn formal_charge(category: &Category, row: usize) -> Result<i8, Diagnostic> {
    if unknown_marker(category, "charge", row) {
        return Ok(0);
    }
    let value = required_integer(category, "charge", row)?;
    i8::try_from(value).map_err(|_| invalid_value(category, "charge", row))
}
