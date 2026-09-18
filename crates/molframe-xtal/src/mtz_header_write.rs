//! Deterministic MTZ header emission.

use super::{
    ReflectionDataset, ReflectionError, ReflectionTable, ReflectionValue, UnitCell,
    calculate_resolution, column_min_max, format_cell, mtz_type, push_record,
};

pub(super) fn write_crystal_header(
    output: &mut Vec<u8>,
    table: &ReflectionTable,
    cell: UnitCell,
    number: i32,
    name: &str,
) -> Result<(), ReflectionError> {
    push_record(output, "VERS MTZ:V1.1");
    push_record(output, &format!("TITLE {}", table.title));
    push_record(
        output,
        &format!("NCOL {} {} 0", table.columns.len(), table.row_count()),
    );
    push_record(output, &format_cell("CELL", None, cell));
    push_record(
        output,
        &format!(
            "SORT {} {} {} {} {}",
            table.sort_order[0],
            table.sort_order[1],
            table.sort_order[2],
            table.sort_order[3],
            table.sort_order[4]
        ),
    );
    let lattice = match name.chars().find(|character| !character.is_whitespace()) {
        Some(value) => value,
        None => 'P',
    };
    push_record(
        output,
        &format!(
            "SYMINF {} {} {} {} '{}' PG?",
            table.symmetry_operations.len(),
            table.symmetry_operations.len(),
            lattice,
            number,
            name
        ),
    );
    for operation in &table.symmetry_operations {
        push_record(output, &format!("SYMM {operation}"));
    }
    let resolution = match table
        .resolution_range
        .or_else(|| calculate_resolution(table, cell).ok())
    {
        Some(value) => value,
        None => [0.0, 0.0],
    };
    push_record(
        output,
        &format!("RESO {:.12} {:.12}", resolution[0], resolution[1]),
    );
    match table.missing_value {
        Some(value) if value.is_finite() => push_record(output, &format!("VALM {value:.9}")),
        Some(_) => return Err(ReflectionError::InvalidNumber),
        None => push_record(output, "VALM NAN"),
    }
    Ok(())
}

pub(super) fn write_column_headers(
    output: &mut Vec<u8>,
    table: &ReflectionTable,
) -> Result<(), ReflectionError> {
    for column in &table.columns {
        if column
            .values
            .iter()
            .any(|value| matches!(value, ReflectionValue::Text(_)))
        {
            return Err(ReflectionError::Unsupported(
                "text-valued MTZ column".into(),
            ));
        }
        let kind = match column.mtz_type {
            Some(value) => value,
            None => mtz_type(column.column_type),
        };
        let (minimum, maximum) = column_min_max(column)?;
        push_record(
            output,
            &format!(
                "COLUMN {:<30} {} {:.9} {:.9} {}",
                column.label, kind, minimum, maximum, column.dataset_id
            ),
        );
    }
    Ok(())
}

pub(super) fn write_dataset_headers(
    output: &mut Vec<u8>,
    table: &ReflectionTable,
    cell: UnitCell,
) -> Result<(), ReflectionError> {
    let datasets = if table.datasets.is_empty() {
        vec![ReflectionDataset {
            id: 0,
            project: "HKL_base".into(),
            crystal: "HKL_base".into(),
            name: "HKL_base".into(),
            wavelength: None,
            cell: Some(cell),
        }]
    } else {
        table.datasets.clone()
    };
    push_record(output, &format!("NDIF {}", datasets.len()));
    for dataset in datasets {
        push_record(
            output,
            &format!("PROJECT {} {}", dataset.id, dataset.project),
        );
        push_record(
            output,
            &format!("CRYSTAL {} {}", dataset.id, dataset.crystal),
        );
        push_record(output, &format!("DATASET {} {}", dataset.id, dataset.name));
        let dataset_cell = match dataset.cell {
            Some(value) => value,
            None => cell,
        };
        push_record(
            output,
            &format_cell("DCELL", Some(dataset.id), dataset_cell),
        );
        if let Some(wavelength) = dataset.wavelength {
            if !wavelength.is_finite() || wavelength <= 0.0 {
                return Err(ReflectionError::InvalidNumber);
            }
            push_record(output, &format!("DWAVEL {} {:.8}", dataset.id, wavelength));
        }
    }
    Ok(())
}

pub(super) fn write_header_tail(output: &mut Vec<u8>, table: &ReflectionTable) {
    for record in &table.extra_header_records {
        push_record(output, record);
    }
    push_record(output, "END");
    if !table.history.is_empty() {
        push_record(output, &format!("MTZHIST {}", table.history.len()));
        for line in &table.history {
            push_record(output, line);
        }
    }
    push_record(output, "MTZENDOFHEADERS");
}
