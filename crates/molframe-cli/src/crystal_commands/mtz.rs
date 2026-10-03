//! Reflection file summaries.

use crate::MtzTable;
use crate::exit::Exit;
use crate::report::{Context, emit_rows};
use molframe::crystal::{ReflectionColumnType, ReflectionTable, read_mtz};
use std::path::Path;

pub(super) fn info(path: &Path, table: MtzTable, context: Context) -> Exit {
    let bytes = match std::fs::read(path) {
        Ok(bytes) => bytes,
        Err(error) => {
            eprintln!("{} could not be read: {error}", path.display());
            return Exit::Failure;
        }
    };
    let reflections = match read_mtz(&bytes) {
        Ok(reflections) => reflections,
        Err(error) => {
            eprintln!("MTZ read failed: {error}");
            return Exit::Consistency;
        }
    };
    match table {
        MtzTable::Summary => emit_rows(
            context,
            &["property", "value"],
            summary(&reflections),
            |(property, value)| vec![property.to_owned(), value],
        ),
        MtzTable::Columns => emit_rows(
            context,
            &["label", "type", "dataset"],
            reflections.columns.iter(),
            |column| {
                vec![
                    column.label.to_string(),
                    column_type(column.column_type).to_owned(),
                    column.dataset_id.to_string(),
                ]
            },
        ),
    }
}

fn summary(table: &ReflectionTable) -> Vec<(&'static str, String)> {
    let mut rows = vec![
        ("title", table.title.to_string()),
        ("reflections", table.row_count().to_string()),
        ("columns", table.columns.len().to_string()),
        ("datasets", table.datasets.len().to_string()),
    ];
    if let Some(number) = table.space_group_number {
        rows.push(("space_group_number", number.to_string()));
    }
    if let Some(name) = &table.space_group_name {
        rows.push(("space_group_name", name.to_string()));
    }
    if let Some(cell) = table.cell {
        for (axis, length) in ["a", "b", "c"].into_iter().zip(cell.lengths) {
            rows.push((cell_property(axis), length.to_string()));
        }
        for (angle, value) in ["alpha", "beta", "gamma"].into_iter().zip(cell.angles) {
            rows.push((cell_property(angle), value.to_string()));
        }
    }
    if let Some([low, high]) = table.resolution_range {
        rows.push(("resolution_low_angstrom", (1.0 / low.sqrt()).to_string()));
        rows.push(("resolution_high_angstrom", (1.0 / high.sqrt()).to_string()));
    }
    rows
}

const fn cell_property(name: &'static str) -> &'static str {
    match name.as_bytes() {
        b"a" => "cell_a",
        b"b" => "cell_b",
        b"c" => "cell_c",
        b"alpha" => "cell_alpha",
        b"beta" => "cell_beta",
        _ => "cell_gamma",
    }
}

const fn column_type(kind: ReflectionColumnType) -> &'static str {
    match kind {
        ReflectionColumnType::MillerIndex => "miller-index",
        ReflectionColumnType::Amplitude => "amplitude",
        ReflectionColumnType::Intensity => "intensity",
        ReflectionColumnType::StandardDeviation => "standard-deviation",
        ReflectionColumnType::Phase => "phase",
        ReflectionColumnType::Flag => "flag",
        ReflectionColumnType::Real => "real",
        ReflectionColumnType::Text => "text",
    }
}
