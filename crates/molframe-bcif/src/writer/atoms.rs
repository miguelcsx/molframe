//! Direct atom-site column encoding without row or document materialisation.

use super::column::{FloatColumnBuilder, IntegerColumnBuilder, TextColumnBuilder};
use super::structure::projection_diagnostic;
use crate::container::EncodedColumn;
use molframe_cif::{CanonicalAtomRow, CanonicalProjection, CanonicalValue};
use molframe_core::diagnostic::Diagnostic;

pub(super) const COLUMN_COUNT: usize = 20;

pub(super) fn visit_columns(
    projection: CanonicalProjection<'_>,
    rows: usize,
    mut visit: impl FnMut(EncodedColumn) -> Result<(), Diagnostic>,
) -> Result<(), Diagnostic> {
    visit(text_column(projection, rows, "group_PDB", |row| {
        CanonicalValue::Present(row.group())
    })?)?;
    visit(integer_column(projection, rows, "id", |row| {
        CanonicalValue::Present(row.atom_site_id())
    })?)?;
    visit(element_column(projection, rows)?)?;
    visit(text_column(projection, rows, "label_atom_id", |row| {
        CanonicalValue::Present(row.label_atom_id())
    })?)?;
    visit(text_column(projection, rows, "label_alt_id", alt_id)?)?;
    visit(text_column(projection, rows, "label_comp_id", |row| {
        CanonicalValue::Present(row.label_comp_id())
    })?)?;
    visit(text_column(projection, rows, "label_asym_id", |row| {
        CanonicalValue::Present(row.label_asym_id())
    })?)?;
    visit(text_column(projection, rows, "label_entity_id", entity_id)?)?;
    visit(integer_column(
        projection,
        rows,
        "label_seq_id",
        label_seq_id,
    )?)?;
    visit(text_column(
        projection,
        rows,
        "pdbx_PDB_ins_code",
        insertion_code,
    )?)?;
    visit(float_column(projection, rows, "Cartn_x", |row| {
        row.coordinate(0)
    })?)?;
    visit(float_column(projection, rows, "Cartn_y", |row| {
        row.coordinate(1)
    })?)?;
    visit(float_column(projection, rows, "Cartn_z", |row| {
        row.coordinate(2)
    })?)?;
    visit(float_column(projection, rows, "occupancy", occupancy)?)?;
    visit(float_column(projection, rows, "B_iso_or_equiv", b_factor)?)?;
    visit(integer_column(
        projection,
        rows,
        "auth_seq_id",
        auth_seq_id,
    )?)?;
    visit(text_column(projection, rows, "auth_comp_id", auth_comp_id)?)?;
    visit(text_column(projection, rows, "auth_asym_id", auth_asym_id)?)?;
    visit(text_column(projection, rows, "auth_atom_id", auth_atom_id)?)?;
    visit(integer_column(
        projection,
        rows,
        "pdbx_PDB_model_num",
        |row| CanonicalValue::Present(row.model_number()),
    )?)
}

// Free accessor functions, not method references: the row methods take `self`
// with the impl's fixed lifetime, so `CanonicalAtomRow::method` fails the
// higher-ranked `for<'row>` bound. A free `fn` late-binds `'row` and satisfies
// it (and reads as a plain function item rather than a redundant closure).
fn alt_id(row: CanonicalAtomRow<'_>) -> CanonicalValue<&str> {
    row.label_alt_id()
}

fn entity_id(row: CanonicalAtomRow<'_>) -> CanonicalValue<&str> {
    row.label_entity_id()
}

fn insertion_code(row: CanonicalAtomRow<'_>) -> CanonicalValue<&str> {
    row.insertion_code()
}

fn auth_comp_id(row: CanonicalAtomRow<'_>) -> CanonicalValue<&str> {
    row.auth_comp_id()
}

fn auth_asym_id(row: CanonicalAtomRow<'_>) -> CanonicalValue<&str> {
    row.auth_asym_id()
}

fn auth_atom_id(row: CanonicalAtomRow<'_>) -> CanonicalValue<&str> {
    row.auth_atom_id()
}

fn label_seq_id(row: CanonicalAtomRow<'_>) -> CanonicalValue<i64> {
    row.label_seq_id()
}

fn auth_seq_id(row: CanonicalAtomRow<'_>) -> CanonicalValue<i64> {
    row.auth_seq_id()
}

fn occupancy(row: CanonicalAtomRow<'_>) -> CanonicalValue<f64> {
    row.occupancy()
}

fn b_factor(row: CanonicalAtomRow<'_>) -> CanonicalValue<f64> {
    row.b_factor()
}

fn text_column<F>(
    projection: CanonicalProjection<'_>,
    rows: usize,
    name: &'static str,
    mut value: F,
) -> Result<EncodedColumn, Diagnostic>
where
    F: for<'row> FnMut(CanonicalAtomRow<'row>) -> CanonicalValue<&'row str>,
{
    let mut column = TextColumnBuilder::new(name, rows);
    projection
        .visit_atom_rows(|row| column.push(value(row)))
        .map_err(|error| projection_diagnostic(&error))?;
    column.finish()
}

fn integer_column<F>(
    projection: CanonicalProjection<'_>,
    rows: usize,
    name: &'static str,
    mut value: F,
) -> Result<EncodedColumn, Diagnostic>
where
    F: for<'row> FnMut(CanonicalAtomRow<'row>) -> CanonicalValue<i64>,
{
    let mut column = IntegerColumnBuilder::new(name, rows);
    projection
        .visit_atom_rows(|row| column.push(value(row)))
        .map_err(|error| projection_diagnostic(&error))?;
    column.finish()
}

fn float_column<F>(
    projection: CanonicalProjection<'_>,
    rows: usize,
    name: &'static str,
    mut value: F,
) -> Result<EncodedColumn, Diagnostic>
where
    F: for<'row> FnMut(CanonicalAtomRow<'row>) -> CanonicalValue<f64>,
{
    let mut column = FloatColumnBuilder::new(name, rows);
    projection
        .visit_atom_rows(|row| column.push(value(row)))
        .map_err(|error| projection_diagnostic(&error))?;
    column.finish()
}

fn element_column(
    projection: CanonicalProjection<'_>,
    rows: usize,
) -> Result<EncodedColumn, Diagnostic> {
    let mut column = TextColumnBuilder::new("type_symbol", rows);
    projection
        .visit_atom_rows(|row| column.push_element(row.element()))
        .map_err(|error| projection_diagnostic(&error))?;
    column.finish()
}
