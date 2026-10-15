//! Comparison output rendering: one row shape per metric family.
//!
//! These are the projections the comparison commands share. They live apart
//! from the commands themselves because they are a single concern — turning a
//! computed value into the requested output kind — and because the command
//! module is at the crate's file-size ceiling.

use crate::report::{Context, Json, Table};

pub(super) fn emit_rmsd(value: f64, atoms: usize, fitted: bool, context: Context) {
    if context.is_json() {
        let mut object = Json::new();
        object
            .number("rmsd", format!("{value:.4}"))
            .number("atoms", atoms)
            .raw("fitted", if fitted { "true" } else { "false" });
        context.result(&object.finish());
    } else if let Some(delimiter) = context.delimiter() {
        let mut table = Table::new(delimiter, &["rmsd", "atoms", "fitted"]);
        let values = [format!("{value:.4}"), atoms.to_string(), fitted.to_string()];
        table.row(values.iter().map(String::as_str));
        context.result(&table.finish());
    } else {
        let how = if fitted {
            "after fitting"
        } else {
            "as they sit"
        };
        context.result(&format!("rmsd {value:.3} over {atoms} atoms, {how}"));
    }
}

pub(super) fn emit_metrics(context: Context, rows: &[(&str, f64)]) {
    if context.is_json() {
        let mut json = Json::new();
        for (name, value) in rows {
            json.number(name, value);
        }
        context.result(&json.finish());
    } else {
        let delimiter = context.table_delimiter();
        let mut table = Table::new(delimiter, &["metric", "value"]);
        for (name, value) in rows {
            let value = value.to_string();
            table.row([*name, &value]);
        }
        context.result(&table.finish());
    }
}

pub(super) fn emit_mappings(context: Context, rows: &[(&str, String, String, f64)]) {
    if context.is_json() {
        let objects = rows
            .iter()
            .map(|(assignment, reference, model, identity)| {
                let mut json = Json::new();
                json.text("assignment", assignment)
                    .text("reference_chain", reference)
                    .text("model_chain", model)
                    .number("identity", identity);
                json.finish()
            })
            .collect::<Vec<_>>();
        context.result(&context.json_records(&objects));
    } else {
        let delimiter = context.table_delimiter();
        let mut table = Table::new(
            delimiter,
            &["assignment", "reference_chain", "model_chain", "identity"],
        );
        for (assignment, reference, model, identity) in rows {
            let identity = identity.to_string();
            table.row([
                *assignment,
                reference.as_str(),
                model.as_str(),
                identity.as_str(),
            ]);
        }
        context.result(&table.finish());
    }
}
