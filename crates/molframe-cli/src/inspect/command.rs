//! CLI rendering over lossless documents and native sequence projection.

use crate::commands::open;
use crate::exit::Exit;
use crate::report::{Context, Json, Table};
use std::fmt::Write as _;
use std::path::Path;

pub(crate) fn assemblies(input: &Path, context: Context) -> Exit {
    use molframe::AssemblyExt as _;
    let structure = match open(input, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let rows = structure
        .assembly_set()
        .into_iter()
        .flat_map(molframe::AssemblySet::assemblies)
        .map(|assembly| {
            vec![
                assembly.id.to_string(),
                match assembly.details.as_deref() {
                    Some(details) => details.to_owned(),
                    None => String::new(),
                },
                assembly
                    .oligomeric
                    .map_or_else(String::new, |value| value.to_string()),
                assembly.generators.len().to_string(),
            ]
        })
        .collect::<Vec<_>>();
    emit_rows(
        context,
        &["assembly", "details", "oligomeric", "generators"],
        &rows,
    );
    Exit::Success
}

pub(crate) fn show(
    input: &Path,
    category_name: &str,
    limit: Option<usize>,
    context: Context,
) -> Exit {
    let document = match molframe::read_document(input) {
        Ok(document) => document,
        Err(findings) => {
            context.findings(&findings, &input.display().to_string());
            return Exit::of(&findings);
        }
    };
    let categories = document
        .blocks()
        .filter_map(|block| block.category(category_name))
        .collect::<Vec<_>>();
    if categories.is_empty() {
        eprintln!("category {category_name:?} is absent");
        return Exit::Usage;
    }
    let header = categories[0].items().collect::<Vec<_>>();
    let rows = categories
        .iter()
        .flat_map(|category| {
            (0..category.row_count()).map(|row| {
                header
                    .iter()
                    .map(|item| value_text(category.value(item, row)))
                    .collect::<Vec<_>>()
            })
        })
        .take(match limit {
            Some(value) => value,
            None => usize::MAX,
        })
        .collect::<Vec<_>>();
    emit_rows(context, &header, &rows);
    Exit::Success
}

pub(crate) fn categories(input: &Path, context: Context) -> Exit {
    let document = match molframe::read_document(input) {
        Ok(document) => document,
        Err(findings) => {
            context.findings(&findings, &input.display().to_string());
            return Exit::of(&findings);
        }
    };
    let rows = document.blocks().flat_map(|block| {
        block.categories().map(move |category| {
            (
                block.name(),
                category.name(),
                category.row_count(),
                category.len(),
            )
        })
    });
    let rows = rows.collect::<Vec<_>>();
    if context.is_json() {
        let objects = rows
            .iter()
            .map(|(block, category, row_count, items)| {
                let mut json = Json::new();
                json.text("block", block)
                    .text("category", category)
                    .number("rows", row_count)
                    .number("items", items);
                json.finish()
            })
            .collect::<Vec<_>>();
        context.result(&context.json_records(&objects));
    } else if let Some(delimiter) = context.delimiter() {
        let mut table = Table::new(delimiter, &["block", "category", "rows", "items"]);
        for (block, category, row_count, items) in rows {
            let values = [
                block.to_owned(),
                category.to_owned(),
                row_count.to_string(),
                items.to_string(),
            ];
            table.row(values.iter().map(String::as_str));
        }
        context.result(&table.finish());
    } else {
        let mut text = String::new();
        for (block, category, row_count, items) in rows {
            let _ = writeln!(
                text,
                "{block}: {category} ({row_count} rows, {items} items)"
            );
        }
        context.result(text.trim_end());
    }
    Exit::Success
}

pub(crate) fn sequence(input: &Path, ccd: &Path, ccd_version: &str, context: Context) -> Exit {
    let structure = match open(input, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let provider = match crate::chemistry::load_ccd(ccd, ccd_version, context) {
        Ok(provider) => provider,
        Err(exit) => return exit,
    };
    let chains =
        match molframe::compare::chain_sequences(&structure, &provider, context.policy.identifiers)
        {
            Ok(chains) => chains,
            Err(finding) => {
                context.findings(&[finding], &input.display().to_string());
                return Exit::Consistency;
            }
        };
    if context.is_json() {
        let objects = chains
            .iter()
            .map(|chain| {
                let mut json = Json::new();
                json.text("chain", &chain.label)
                    .text("sequence", &String::from_utf8_lossy(&chain.sequence))
                    .number("length", chain.sequence.len());
                json.finish()
            })
            .collect::<Vec<_>>();
        context.result(&context.json_records(&objects));
    } else if let Some(delimiter) = context.delimiter() {
        let mut table = Table::new(delimiter, &["chain", "sequence", "length"]);
        for chain in chains {
            let values = [
                chain.label,
                String::from_utf8_lossy(&chain.sequence).into_owned(),
                chain.sequence.len().to_string(),
            ];
            table.row(values.iter().map(String::as_str));
        }
        context.result(&table.finish());
    } else {
        let mut fasta = String::new();
        for chain in chains {
            let _ = writeln!(fasta, ">{}", chain.label);
            let _ = writeln!(fasta, "{}", String::from_utf8_lossy(&chain.sequence));
        }
        context.result(fasta.trim_end());
    }
    Exit::Success
}

pub(super) fn value_text(value: Option<&molframe::CifValue>) -> String {
    match value {
        Some(molframe::CifValue::Inapplicable) => ".".to_owned(),
        Some(molframe::CifValue::Unknown) | None => "?".to_owned(),
        Some(molframe::CifValue::Text(value)) => value.to_string(),
        Some(molframe::CifValue::Integer(value)) => value.to_string(),
        Some(molframe::CifValue::Float(value)) => value.to_string(),
    }
}

fn emit_rows(context: Context, header: &[&str], rows: &[Vec<String>]) {
    if context.is_json() {
        let objects = rows
            .iter()
            .map(|row| {
                let mut json = Json::new();
                for (key, value) in header.iter().zip(row) {
                    json.text(key, value);
                }
                json.finish()
            })
            .collect::<Vec<_>>();
        context.result(&context.json_records(&objects));
    } else if let Some(delimiter) = context.delimiter() {
        let mut table = Table::new(delimiter, header);
        for row in rows {
            table.row(row.iter().map(String::as_str));
        }
        context.result(&table.finish());
    } else {
        let mut table = Table::new('\t', header);
        for row in rows {
            table.row(row.iter().map(String::as_str));
        }
        context.result(&table.finish());
    }
}
