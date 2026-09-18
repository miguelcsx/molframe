//! Declarative selection evaluation and native structure materialisation.

use crate::commands::open;
use crate::exit::Exit;
use crate::report::{Context, Json, OutputKind, Table};
use molframe::{AtomIndex, ChainIndex, QueryStructure as _};
use std::path::Path;

pub(crate) fn select(
    input: &Path,
    query: &str,
    output: Option<&Path>,
    count_only: bool,
    context: Context,
) -> Exit {
    let structure = match open(input, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let evaluation = match structure.select_text(query, context.policy, context.execution) {
        Ok(evaluation) => evaluation,
        Err(findings) => {
            context.findings(&findings, &input.display().to_string());
            return Exit::of(&findings);
        }
    };
    context.findings(&evaluation.warnings, &input.display().to_string());
    if count_only {
        emit_count(evaluation.selection.len(), context);
        return Exit::Success;
    }
    let selected = match structure.materialize(&evaluation.selection) {
        Ok(selected) => selected,
        Err(findings) => {
            context.findings(&findings, &input.display().to_string());
            return Exit::of(&findings);
        }
    };
    if let Some(output) = output {
        write_selection(output, &selected, context)
    } else {
        emit_atoms(&selected, context);
        Exit::Success
    }
}

fn write_selection(output: &Path, selected: &molframe::Structure, context: Context) -> Exit {
    let table_result = match context.format {
        OutputKind::Arrow => {
            molframe::write_atom_ipc_with_metadata(output, selected, context.provenance_metadata())
        }
        OutputKind::Parquet => molframe::write_atom_parquet_with_metadata(
            output,
            selected,
            context.provenance_metadata(),
        ),
        OutputKind::Text
        | OutputKind::Json
        | OutputKind::JsonLines
        | OutputKind::Csv
        | OutputKind::Tsv => return write_structure(output, selected, context),
    };
    match table_result {
        Ok(()) => {
            emit_written(output, selected, context);
            Exit::Success
        }
        Err(error) => {
            eprintln!("selection export failed: {error}");
            Exit::Failure
        }
    }
}

fn write_structure(output: &Path, selected: &molframe::Structure, context: Context) -> Exit {
    match molframe::write(output, selected) {
        Ok(()) => {
            emit_written(output, selected, context);
            Exit::Success
        }
        Err(findings) => {
            context.findings(&findings, &output.display().to_string());
            Exit::of(&findings)
        }
    }
}

fn emit_written(output: &Path, selected: &molframe::Structure, context: Context) {
    if context.is_json() {
        let mut object = Json::new();
        object
            .text("output", &output.display().to_string())
            .number("atoms", selected.atom_count())
            .number("residues", selected.residue_count())
            .number("chains", selected.chain_count());
        context.result(&object.finish());
    } else {
        context.result(&format!(
            "wrote {} atoms to {}",
            selected.atom_count(),
            output.display()
        ));
    }
}

fn emit_count(count: u64, context: Context) {
    if context.is_json() {
        let mut object = Json::new();
        object.number("atoms", count);
        context.result(&object.finish());
    } else if let Some(delimiter) = context.delimiter() {
        let mut table = Table::new(delimiter, &["atoms"]);
        let count = count.to_string();
        table.row([count.as_str()]);
        context.result(&table.finish());
    } else {
        context.result(&count.to_string());
    }
}

fn emit_atoms(structure: &molframe::Structure, context: Context) {
    let rows: Vec<AtomRow> = (0..structure.atom_count())
        .filter_map(|position| atom_row(structure, AtomIndex::new(position)))
        .collect();
    if context.is_json() {
        let records: Vec<String> = rows.iter().map(AtomRow::json).collect();
        context.result(&context.json_records(&records));
    } else {
        let delimiter = match context.delimiter() {
            Some(delimiter) => delimiter,
            None => '\t',
        };
        let mut table = Table::new(
            delimiter,
            &["atom", "name", "element", "chain", "residue", "x", "y", "z"],
        );
        for row in &rows {
            let values = row.values();
            table.row(values.iter().map(String::as_str));
        }
        context.result(&table.finish());
    }
}

#[derive(Debug)]
struct AtomRow {
    atom: u32,
    name: String,
    element: String,
    chain: String,
    residue: String,
    position: Option<[f32; 3]>,
}

impl AtomRow {
    fn values(&self) -> [String; 8] {
        let coordinates = match self.position {
            Some(position) => position.map(|value| value.to_string()),
            None => std::array::from_fn(|_| String::new()),
        };
        [
            self.atom.to_string(),
            self.name.clone(),
            self.element.clone(),
            self.chain.clone(),
            self.residue.clone(),
            coordinates[0].clone(),
            coordinates[1].clone(),
            coordinates[2].clone(),
        ]
    }

    fn json(&self) -> String {
        let mut object = Json::new();
        object
            .number("atom", self.atom)
            .text("name", &self.name)
            .text("element", &self.element)
            .text("chain", &self.chain)
            .text("residue", &self.residue);
        match self.position {
            Some(position) => {
                object
                    .number("x", position[0])
                    .number("y", position[1])
                    .number("z", position[2]);
            }
            None => {
                object.raw("x", "null").raw("y", "null").raw("z", "null");
            }
        }
        object.finish()
    }
}

fn atom_row(structure: &molframe::Structure, index: AtomIndex) -> Option<AtomRow> {
    let atom = structure.atom(index)?;
    let residue = atom.residue()?;
    let chain_index = structure
        .data()
        .topology
        .chains
        .containing(residue.index().get())?;
    let chain = structure.chain(ChainIndex::new(chain_index.get()))?;
    Some(AtomRow {
        atom: index.get(),
        name: atom.name().map_or_else(String::new, str::to_owned),
        element: atom
            .element()
            .map_or_else(String::new, |element| element.to_string()),
        chain: chain.label().map_or_else(String::new, str::to_owned),
        residue: residue.name().map_or_else(String::new, str::to_owned),
        position: atom.position(),
    })
}
