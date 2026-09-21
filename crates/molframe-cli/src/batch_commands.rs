//! Deterministic bounded parallel execution over path and glob inputs.

use crate::BatchCommand;
use crate::exit::Exit;
use crate::report::{Context, Json, Table};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

pub(crate) fn execute(command: BatchCommand, context: Context) -> Exit {
    match command {
        BatchCommand::Info { inputs } => info(&inputs, context),
        BatchCommand::Convert { inputs, to, outdir } => convert(&inputs, &to, &outdir, context),
    }
}

fn info(patterns: &[String], context: Context) -> Exit {
    let paths = match expand(patterns) {
        Ok(paths) => paths,
        Err(error) => return batch_error(&error),
    };
    let results = run(&paths, context.execution, |path| {
        molframe::read_with_options(
            path,
            &molframe::ReadOptions::new()
                .mode(context.mode)
                .missing_element_policy(context.missing_element_policy)
                .ambiguous_residue_boundary_policy(context.residue_boundary_policy),
        )
        .map(|(structure, _findings)| {
            vec![
                path.display().to_string(),
                structure.model_count().to_string(),
                structure.chain_count().to_string(),
                structure.residue_count().to_string(),
                structure.atom_count().to_string(),
            ]
        })
        .map_err(|findings| findings_message(&findings))
    });
    emit_batch(
        &results,
        &["input", "models", "chains", "residues", "atoms"],
        context,
    )
}

fn convert(patterns: &[String], format: &str, outdir: &Path, context: Context) -> Exit {
    if !matches!(format, "cif" | "bcif" | "pdb") {
        eprintln!("batch convert --to must be cif, bcif or pdb");
        return Exit::Usage;
    }
    let paths = match expand(patterns) {
        Ok(paths) => paths,
        Err(error) => return batch_error(&error),
    };
    if let Err(error) = std::fs::create_dir_all(outdir) {
        eprintln!("could not create {}: {error}", outdir.display());
        return Exit::Failure;
    }
    let destinations = match destinations(&paths, outdir, format) {
        Ok(destinations) => destinations,
        Err(error) => return batch_error(&error),
    };
    let jobs = paths.into_iter().zip(destinations).collect::<Vec<_>>();
    let results = run(&jobs, context.execution, |(input, output)| {
        if output.exists() {
            return Err(format!("destination exists: {}", output.display()));
        }
        let (structure, _findings) = molframe::read_with_options(
            input,
            &molframe::ReadOptions::new()
                .mode(context.mode)
                .missing_element_policy(context.missing_element_policy)
                .ambiguous_residue_boundary_policy(context.residue_boundary_policy),
        )
        .map_err(|findings| findings_message(&findings))?;
        molframe::write(output, &structure).map_err(|findings| findings_message(&findings))?;
        Ok(vec![
            input.display().to_string(),
            output.display().to_string(),
            structure.atom_count().to_string(),
        ])
    });
    emit_batch(&results, &["input", "output", "atoms"], context)
}

fn run<T: Sync, R: Send>(
    jobs: &[T],
    context: &molframe_core::ExecutionContext,
    operation: impl Fn(&T) -> Result<R, String> + Sync + Send,
) -> Vec<Result<R, String>> {
    let plan = molframe_core::parallel::BlockPlan::new(jobs.len(), 1);
    match molframe_core::parallel::map_blocks_in(plan, context, |_, range| {
        let Some(job) = jobs.get(range.start) else {
            return Err("batch plan produced an invalid job index".to_owned());
        };
        operation(job)
    }) {
        Ok(results) => results,
        Err(error) => vec![Err(format!("shared worker pool failed: {error}"))],
    }
}

fn expand(patterns: &[String]) -> Result<Vec<PathBuf>, String> {
    let mut paths = BTreeSet::new();
    for pattern in patterns {
        let entries = glob::glob(pattern).map_err(|error| error.to_string())?;
        for entry in entries {
            let path = entry.map_err(|error| error.to_string())?;
            if path.is_file() {
                paths.insert(path);
            }
        }
    }
    if paths.is_empty() {
        Err("no input files matched".to_owned())
    } else {
        Ok(paths.into_iter().collect())
    }
}

fn destinations(paths: &[PathBuf], outdir: &Path, format: &str) -> Result<Vec<PathBuf>, String> {
    let mut seen = BTreeSet::new();
    let mut output = Vec::with_capacity(paths.len());
    for path in paths {
        let stem = path
            .file_stem()
            .ok_or_else(|| format!("input has no file stem: {}", path.display()))?;
        let destination = outdir.join(stem).with_extension(format);
        if !seen.insert(destination.clone()) {
            return Err(format!("multiple inputs map to {}", destination.display()));
        }
        output.push(destination);
    }
    Ok(output)
}

fn emit_batch(results: &[Result<Vec<String>, String>], header: &[&str], context: Context) -> Exit {
    let mut rows = Vec::new();
    let mut failures = Vec::new();
    for (index, result) in results.iter().enumerate() {
        eprintln!("batch {}/{}", index + 1, results.len());
        match result {
            Ok(row) => rows.push(row.clone()),
            Err(error) => failures.push(error.clone()),
        }
    }
    if context.is_json() {
        let records = rows
            .iter()
            .map(|row| {
                let mut json = Json::new();
                for (name, value) in header.iter().zip(row) {
                    json.text(name, value);
                }
                json.finish()
            })
            .collect::<Vec<_>>();
        let mut envelope = Json::new();
        envelope
            .raw("results", &format!("[{}]", records.join(",")))
            .raw(
                "errors",
                &format!(
                    "[{}]",
                    failures
                        .iter()
                        .map(|error| format!("\"{}\"", error.replace('"', "\\\"")))
                        .collect::<Vec<_>>()
                        .join(",")
                ),
            );
        context.result(&envelope.finish());
    } else {
        let delimiter = context.table_delimiter();
        let mut table = Table::new(delimiter, header);
        for row in &rows {
            table.row(row.iter().map(String::as_str));
        }
        context.result(&table.finish());
        for failure in &failures {
            eprintln!("batch item failed: {failure}");
        }
    }
    if failures.is_empty() {
        Exit::Success
    } else {
        Exit::Failure
    }
}

fn findings_message(findings: &[molframe::Diagnostic]) -> String {
    findings
        .iter()
        .map(|finding| finding.code().to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn batch_error(message: &str) -> Exit {
    eprintln!("batch refused: {message}");
    Exit::Usage
}
