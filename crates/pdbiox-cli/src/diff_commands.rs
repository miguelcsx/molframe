//! Rendering for the native semantic structure comparison.

use crate::commands::open;
use crate::exit::Exit;
use crate::report::{Context, Json, Table};
use pdbiox::{CountDifference, StructureDifference, ValueDifference};
use std::fmt::Write as _;
use std::path::Path;

pub(crate) fn diff(
    left_path: &Path,
    right_path: &Path,
    coordinate_tolerance: f32,
    context: Context,
) -> Exit {
    let left = match open(left_path, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let right = match open(right_path, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let result = pdbiox::structure_difference(
        &left,
        &right,
        pdbiox::StructureDifferenceOptions {
            coordinate_tolerance,
        },
    );
    let difference = match result {
        Ok(difference) => difference,
        Err(error) => {
            eprintln!("comparison refused: {error}");
            return Exit::Usage;
        }
    };

    if context.is_json() {
        emit_json(&difference, coordinate_tolerance, context);
    } else if let Some(delimiter) = context.delimiter() {
        context.result(&emit_table(&difference, delimiter));
    } else {
        context.result(&emit_text(&difference, coordinate_tolerance));
    }
    Exit::Success
}

fn emit_json(difference: &StructureDifference, tolerance: f32, context: Context) {
    let mut object = Json::new();
    object
        .raw("equivalent", bool_json(difference.is_empty()))
        .number("coordinate_tolerance_angstrom", tolerance)
        .number("changed_entities", difference.changed_entities)
        .number("changed_chains", difference.changed_chains)
        .number("changed_residues", difference.changed_residues)
        .number("changed_atoms", difference.changed_atoms)
        .number("changed_bonds", difference.changed_bonds)
        .number("changed_positions", difference.changed_positions);
    if let Some(maximum) = difference.maximum_displacement {
        object.number("maximum_displacement_angstrom", maximum);
    }
    add_count_json(&mut object, "models", difference.models.as_ref());
    add_count_json(&mut object, "entities", difference.entities.as_ref());
    add_count_json(&mut object, "chains", difference.chains.as_ref());
    add_count_json(&mut object, "residues", difference.residues.as_ref());
    add_count_json(&mut object, "atoms", difference.atoms.as_ref());
    add_count_json(&mut object, "bonds", difference.bonds.as_ref());
    object.raw("metadata", &metadata_json(difference));
    context.result(&object.finish());
}

fn add_count_json(object: &mut Json, name: &str, difference: Option<&CountDifference>) {
    if let Some(difference) = difference {
        object.raw(
            &format!("{name}_count"),
            &format!(
                "{{\"left\":{},\"right\":{}}}",
                difference.left, difference.right
            ),
        );
    }
}

fn metadata_json(difference: &StructureDifference) -> String {
    let mut fields = Vec::new();
    push_optional_text_json(&mut fields, "id", difference.metadata.id.as_ref());
    push_optional_text_json(&mut fields, "title", difference.metadata.title.as_ref());
    push_optional_text_json(&mut fields, "method", difference.metadata.method.as_ref());
    push_optional_number_json(
        &mut fields,
        "resolution",
        difference.metadata.resolution.as_ref(),
    );
    push_cell_json(&mut fields, difference.metadata.cell.as_ref());
    format!("[{}]", fields.join(","))
}

fn push_optional_text_json(
    output: &mut Vec<String>,
    field: &str,
    difference: Option<&ValueDifference<Option<Box<str>>>>,
) {
    if let Some(difference) = difference {
        let mut object = Json::new();
        object.text("field", field);
        optional_text_field(&mut object, "left", difference.left.as_deref());
        optional_text_field(&mut object, "right", difference.right.as_deref());
        output.push(object.finish());
    }
}

fn optional_text_field(object: &mut Json, name: &str, value: Option<&str>) {
    match value {
        Some(value) => {
            object.text(name, value);
        }
        None => {
            object.raw(name, "null");
        }
    }
}

fn push_optional_number_json(
    output: &mut Vec<String>,
    field: &str,
    difference: Option<&ValueDifference<Option<f32>>>,
) {
    if let Some(difference) = difference {
        let mut object = Json::new();
        object.text("field", field);
        optional_number_field(&mut object, "left", difference.left);
        optional_number_field(&mut object, "right", difference.right);
        output.push(object.finish());
    }
}

fn optional_number_field(object: &mut Json, name: &str, value: Option<f32>) {
    match value {
        Some(value) => {
            object.number(name, value);
        }
        None => {
            object.raw(name, "null");
        }
    }
}

fn push_cell_json(
    output: &mut Vec<String>,
    difference: Option<&ValueDifference<Option<pdbiox::UnitCell>>>,
) {
    if let Some(difference) = difference {
        let mut object = Json::new();
        object.text("field", "cell");
        cell_field(&mut object, "left", difference.left);
        cell_field(&mut object, "right", difference.right);
        output.push(object.finish());
    }
}

fn cell_field(object: &mut Json, name: &str, cell: Option<pdbiox::UnitCell>) {
    match cell {
        Some(cell) => object.raw(
            name,
            &format!(
                "{{\"lengths\":[{},{},{}],\"angles\":[{},{},{}]}}",
                cell.lengths[0],
                cell.lengths[1],
                cell.lengths[2],
                cell.angles[0],
                cell.angles[1],
                cell.angles[2]
            ),
        ),
        None => object.raw(name, "null"),
    };
}

fn emit_table(difference: &StructureDifference, delimiter: char) -> String {
    let mut table = Table::new(delimiter, &["scope", "field", "left", "right"]);
    push_count_rows(&mut table, difference);
    push_metadata_rows(&mut table, difference);
    let summary = [
        "summary".to_owned(),
        "changed_positions".to_owned(),
        difference.changed_positions.to_string(),
        String::new(),
    ];
    table.row(summary.iter().map(String::as_str));
    table.finish()
}

fn push_count_rows(table: &mut Table, difference: &StructureDifference) {
    for (field, value) in [
        ("models", difference.models.as_ref()),
        ("entities", difference.entities.as_ref()),
        ("chains", difference.chains.as_ref()),
        ("residues", difference.residues.as_ref()),
        ("atoms", difference.atoms.as_ref()),
        ("bonds", difference.bonds.as_ref()),
    ] {
        if let Some(value) = value {
            let row = [
                "count".to_owned(),
                field.to_owned(),
                value.left.to_string(),
                value.right.to_string(),
            ];
            table.row(row.iter().map(String::as_str));
        }
    }
}

fn push_metadata_rows(table: &mut Table, difference: &StructureDifference) {
    push_value_row(table, "id", difference.metadata.id.as_ref());
    push_value_row(table, "title", difference.metadata.title.as_ref());
    push_value_row(table, "method", difference.metadata.method.as_ref());
    push_value_row(table, "resolution", difference.metadata.resolution.as_ref());
    push_value_row(table, "cell", difference.metadata.cell.as_ref());
}

fn push_value_row<T: std::fmt::Debug>(
    table: &mut Table,
    field: &str,
    difference: Option<&ValueDifference<T>>,
) {
    if let Some(difference) = difference {
        let row = [
            "metadata".to_owned(),
            field.to_owned(),
            format!("{:?}", difference.left),
            format!("{:?}", difference.right),
        ];
        table.row(row.iter().map(String::as_str));
    }
}

fn emit_text(difference: &StructureDifference, tolerance: f32) -> String {
    if difference.is_empty() {
        return format!("equivalent (coordinate tolerance {tolerance} angstrom)");
    }
    let mut text = String::from("structures differ");
    write_count_text(&mut text, "models", difference.models.as_ref());
    write_count_text(&mut text, "entities", difference.entities.as_ref());
    write_count_text(&mut text, "chains", difference.chains.as_ref());
    write_count_text(&mut text, "residues", difference.residues.as_ref());
    write_count_text(&mut text, "atoms", difference.atoms.as_ref());
    write_count_text(&mut text, "bonds", difference.bonds.as_ref());
    let _ = write!(
        text,
        "\nchanged rows: entities {}, chains {}, residues {}, atoms {}, bonds {}",
        difference.changed_entities,
        difference.changed_chains,
        difference.changed_residues,
        difference.changed_atoms,
        difference.changed_bonds
    );
    let _ = write!(
        text,
        "\nchanged positions: {}",
        difference.changed_positions
    );
    if let Some(maximum) = difference.maximum_displacement {
        let _ = write!(text, " (maximum {maximum} angstrom)");
    }
    let metadata = metadata_names(difference);
    if !metadata.is_empty() {
        let _ = write!(text, "\nchanged metadata: {}", metadata.join(", "));
    }
    text
}

fn write_count_text(output: &mut String, name: &str, difference: Option<&CountDifference>) {
    if let Some(difference) = difference {
        let _ = write!(
            output,
            "\n{name}: {} -> {}",
            difference.left, difference.right
        );
    }
}

fn metadata_names(difference: &StructureDifference) -> Vec<&'static str> {
    [
        ("id", difference.metadata.id.is_some()),
        ("title", difference.metadata.title.is_some()),
        ("method", difference.metadata.method.is_some()),
        ("resolution", difference.metadata.resolution.is_some()),
        ("cell", difference.metadata.cell.is_some()),
    ]
    .into_iter()
    .filter_map(|(name, changed)| changed.then_some(name))
    .collect()
}

const fn bool_json(value: bool) -> &'static str {
    if value { "true" } else { "false" }
}
