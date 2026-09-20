//! What each subcommand does.
//!
//! Every one of these reads through the library, formats what came back, and
//! returns an exit code. None of them computes anything: a convenience worth
//! having belongs in the library, where the Rust callers get it too.

use crate::exit::Exit;
use crate::report::{Context, Json, Table, json_array};
use molframe::formats::pdb::{PdbIdentifierNamespace, PdbOptions};
use molframe::query::{Evaluation, Groups};
use molframe::{ChainRef, Findings, Format, Namespace, Query, ReadOptions, Structure};
use std::fmt::Write as _;
use std::path::Path;

/// Compiles textual syntax and evaluates it under the CLI's execution context.
///
/// The CLI governs resources explicitly (memory, workers, cancellation), so
/// selection goes through the full-control verb rather than the defaulting
/// `QueryStructure::select_text`.
pub(super) fn select_text(
    structure: &Structure,
    source: &str,
    policy: &molframe::AnalysisPolicy,
    execution: &molframe_core::ExecutionContext,
) -> Result<Evaluation, Findings> {
    use molframe::QueryStructure as _;
    let query = Query::compile(source).map_err(Findings::from)?;
    structure.select_with_options(&query, policy, &Groups::new(), execution)
}

/// Reads a file, printing whatever was wrong with it.
pub(super) fn open(path: &Path, context: Context) -> Result<Structure, Exit> {
    let options = ReadOptions::new()
        .mode(context.mode)
        .missing_element_policy(context.missing_element_policy)
        .ambiguous_residue_boundary_policy(context.residue_boundary_policy);
    let origin = path.display().to_string();
    match molframe::read_with_options(path, &options) {
        Ok((structure, findings)) => {
            context.findings(&findings, &origin);
            Ok(structure)
        }
        Err(findings) => {
            context.findings(&findings, &origin);
            Err(Exit::of(&findings))
        }
    }
}

/// Summarises a structure.
pub fn info(path: &Path, detail: bool, context: Context) -> Exit {
    let structure = match open(path, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let namespace = match context.policy.identifiers {
        Namespace::Label => Namespace::Label,
        Namespace::Auth => Namespace::Auth,
        Namespace::Explicit => {
            eprintln!(
                "info requires label or auth identifiers; explicit mappings are not available"
            );
            return Exit::Usage;
        }
        _ => {
            eprintln!("info does not support the requested identifier namespace");
            return Exit::Usage;
        }
    };
    let data = structure.engine().data();

    if context.is_json() {
        let mut object = Json::new();
        if let Some(id) = &data.entry.id {
            object.text("id", id);
        }
        object
            .number("models", structure.model_count())
            .number("chains", structure.chain_count())
            .number("entities", structure.engine().entity_count())
            .number("residues", structure.residue_count())
            .number("atoms", structure.atom_count());
        if detail {
            object.raw("chain_detail", &chain_detail_json(&structure, namespace));
        }
        context.result(&object.finish());
        return Exit::Success;
    }
    if let Some(delimiter) = context.delimiter() {
        context.result(&info_table(&structure, detail, delimiter, namespace));
        return Exit::Success;
    }

    let mut text = String::new();
    if let Some(id) = &data.entry.id {
        let _ = writeln!(text, "entry     {id}");
    }
    if let Some(title) = &data.entry.title {
        let _ = writeln!(text, "title     {title}");
    }
    let _ = writeln!(text, "models    {}", structure.model_count());
    let _ = writeln!(text, "chains    {}", structure.chain_count());
    let _ = writeln!(text, "entities  {}", structure.engine().entity_count());
    let _ = writeln!(text, "residues  {}", structure.residue_count());
    let _ = write!(text, "atoms     {}", structure.atom_count());
    if let Some(cell) = data.cell {
        let _ = write!(
            text,
            "\ncell      {:.3} {:.3} {:.3}  {:.2} {:.2} {:.2}",
            cell.lengths[0],
            cell.lengths[1],
            cell.lengths[2],
            cell.angles[0],
            cell.angles[1],
            cell.angles[2],
        );
    }
    if detail {
        for chain in data.chains() {
            let _ = write!(
                text,
                "\n  chain {:<4} residues {:>6}  atoms {:>7}",
                label_of(&structure, chain, namespace),
                chain.residues().count(),
                chain
                    .residues()
                    .map(|residue| residue.atoms().count())
                    .sum::<usize>(),
            );
        }
    }
    context.result(&text);
    Exit::Success
}

fn info_table(
    structure: &Structure,
    detail: bool,
    delimiter: char,
    namespace: Namespace,
) -> String {
    let id = match &structure.engine().data().entry.id {
        Some(id) => id.as_ref(),
        None => "",
    };
    if detail {
        let mut table = Table::new(delimiter, &["entry", "chain", "residues", "atoms"]);
        for chain in structure.engine().data().chains() {
            let values = [
                id.to_owned(),
                label_of(structure, chain, namespace),
                chain.residues().count().to_string(),
                chain
                    .residues()
                    .map(|residue| residue.atoms().count())
                    .sum::<usize>()
                    .to_string(),
            ];
            table.row(values.iter().map(String::as_str));
        }
        return table.finish();
    }
    let mut table = Table::new(
        delimiter,
        &["entry", "models", "chains", "entities", "residues", "atoms"],
    );
    let values = [
        id.to_owned(),
        structure.model_count().to_string(),
        structure.chain_count().to_string(),
        structure.engine().entity_count().to_string(),
        structure.residue_count().to_string(),
        structure.atom_count().to_string(),
    ];
    table.row(values.iter().map(String::as_str));
    table.finish()
}

fn chain_detail_json(structure: &Structure, namespace: Namespace) -> String {
    let entries: Vec<String> = structure
        .engine()
        .data()
        .chains()
        .map(|chain| {
            let mut object = Json::new();
            object
                .text("chain", &label_of(structure, chain, namespace))
                .number("residues", chain.residues().count())
                .number(
                    "atoms",
                    chain
                        .residues()
                        .map(|residue| residue.atoms().count())
                        .sum::<usize>(),
                );
            object.finish()
        })
        .collect();
    json_array(&entries)
}

fn label_of(structure: &Structure, chain: ChainRef<'_>, namespace: Namespace) -> String {
    let symbol = match namespace {
        Namespace::Label => chain.label_asym_id(),
        Namespace::Auth => chain.auth_asym_id(),
        _ => None,
    };
    match symbol.and_then(|symbol| structure.engine().resolve(symbol)) {
        Some(label) => label.to_owned(),
        None => "?".to_owned(),
    }
}

/// Converts a structure from one format to another.
pub fn convert(
    input: &Path,
    output: &Path,
    chain_map: &[String],
    hybrid36: bool,
    preserve: bool,
    cif_options: &molframe::formats::cif::CifWriteOptions,
    context: Context,
) -> Exit {
    if preserve {
        return preserving_convert(input, output, context);
    }
    let structure = match open(input, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };

    // The output's name decides the format, because that is what the caller
    // wrote down; nothing is inferred from the input.
    let Some(target) = output
        .file_name()
        .and_then(|name| name.to_str())
        .and_then(Format::from_name)
    else {
        eprintln!(
            "cannot tell what format {} should be; use a recognised extension",
            output.display()
        );
        return Exit::Usage;
    };

    let pdb_options = if matches!(target, Format::Pdb | Format::Pqr | Format::Pdbqt) {
        match pdb_options(chain_map, hybrid36, context.policy.identifiers) {
            Ok(options) => Some(options),
            Err(exit) => return exit,
        }
    } else {
        None
    };
    if !matches!(
        target,
        Format::Mmcif
            | Format::BinaryCif
            | Format::Mmtf
            | Format::Pdb
            | Format::Pqr
            | Format::Pdbqt
    ) {
        return Exit::Usage;
    }
    let mut options = molframe::WriteOptions::canonical().with_cif(cif_options.clone());
    if let Some(pdb) = pdb_options.clone() {
        options = options.with_pdb(pdb);
    }
    let written = molframe::write_with_options(output, &structure, &options);
    match written {
        Ok(()) => Exit::Success,
        Err(findings) => {
            context.findings(&findings, &input.display().to_string());
            Exit::of(&findings)
        }
    }
}

fn pdb_options(
    chain_map: &[String],
    hybrid36: bool,
    namespace: Namespace,
) -> Result<PdbOptions, Exit> {
    let namespace = match namespace {
        Namespace::Label => PdbIdentifierNamespace::Label,
        Namespace::Auth => PdbIdentifierNamespace::Auth,
        Namespace::Explicit => {
            eprintln!("PDB output requires label or auth identifiers, not explicit namespace");
            return Err(Exit::Usage);
        }
        _ => {
            eprintln!("PDB output does not support the requested identifier namespace");
            return Err(Exit::Usage);
        }
    };
    let mut options = PdbOptions::new().hybrid36(hybrid36).namespace(namespace);
    for mapping in chain_map.iter().flat_map(|value| value.split(',')) {
        let Some((from, to)) = mapping.split_once('=') else {
            eprintln!("a chain mapping must be written FROM=TO, not {mapping:?}");
            return Err(Exit::Usage);
        };
        options = options.chain_map(from, to);
    }
    Ok(options)
}

fn preserving_convert(input: &Path, output: &Path, context: Context) -> Exit {
    if output
        .file_name()
        .and_then(|name| name.to_str())
        .and_then(Format::from_name)
        != Some(Format::Mmcif)
    {
        eprintln!("--preserve requires an mmCIF destination");
        return Exit::Usage;
    }
    let document = match molframe::read_document(input) {
        Ok(document) => document,
        Err(findings) => {
            context.findings(&findings, &input.display().to_string());
            return Exit::of(&findings);
        }
    };
    let mut sink =
        match molframe_core::io::OutputSink::create(output, molframe::OutputOptions::default()) {
            Ok(sink) => sink,
            Err(finding) => {
                context.findings(&[finding], &output.display().to_string());
                return Exit::Failure;
            }
        };
    if let Err(error) = molframe::formats::cif::write_preserving_to(&document, &mut sink) {
        eprintln!("could not write {}: {error}", output.display());
        return Exit::Failure;
    }
    match sink.finish() {
        Ok(()) => Exit::Success,
        Err(finding) => {
            context.findings(&[finding], &output.display().to_string());
            Exit::Failure
        }
    }
}

/// Measures a structure's overall shape.
pub fn measure(path: &Path, context: Context) -> Exit {
    let structure = match open(path, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let positions = structure.coordinates();

    let Some(centre) = molframe::geometry::centroid(positions) else {
        eprintln!("the structure holds no positions to measure");
        return Exit::Consistency;
    };
    let radius = molframe::geometry::radius_of_gyration(positions, &[]);

    if context.is_json() {
        let mut object = Json::new();
        object.number("atoms", positions.len());
        for (axis, value) in ["x", "y", "z"].into_iter().zip(centre) {
            object.number(&format!("centre_{axis}"), format!("{value:.4}"));
        }
        if let Some(radius) = radius {
            object.number("radius_of_gyration", format!("{radius:.4}"));
        }
        context.result(&object.finish());
    } else if let Some(delimiter) = context.delimiter() {
        let values = [
            positions.len().to_string(),
            format!("{:.4}", centre[0]),
            format!("{:.4}", centre[1]),
            format!("{:.4}", centre[2]),
            match radius {
                Some(radius) => format!("{radius:.4}"),
                None => String::new(),
            },
        ];
        context.result(&one_row(
            delimiter,
            &[
                "atoms",
                "centre_x",
                "centre_y",
                "centre_z",
                "radius_of_gyration",
            ],
            &values,
        ));
    } else {
        let mut text = format!("atoms     {}\n", positions.len());
        let _ = write!(
            text,
            "centre    {:.3} {:.3} {:.3}",
            centre[0], centre[1], centre[2]
        );
        if let Some(radius) = radius {
            let _ = write!(text, "\nradius    {radius:.3}");
        }
        context.result(&text);
    }
    Exit::Success
}

/// Prints the effective policy after configuration and command-line overrides.
pub fn policy(context: Context) -> Exit {
    let policy = context.policy;
    if context.is_json() {
        let mut object = Json::new();
        object.text(
            "profile",
            &match policy.profile() {
                Some(profile) => profile.to_string(),
                None => "(modified)".to_owned(),
            },
        );
        object.text("fingerprint", &policy.fingerprint().to_string());
        context.result(&object.finish());
    } else if let Some(delimiter) = context.delimiter() {
        let values = [
            match policy.profile() {
                Some(profile) => profile.to_string(),
                None => "(modified)".to_owned(),
            },
            policy.fingerprint().to_string(),
        ];
        context.result(&one_row(delimiter, &["profile", "fingerprint"], &values));
    } else {
        context.result(&format!(
            "{policy}\n\nfingerprint      {}",
            policy.fingerprint()
        ));
    }
    Exit::Success
}

fn one_row(delimiter: char, header: &[&str], values: &[String]) -> String {
    let mut table = Table::new(delimiter, header);
    table.row(values.iter().map(String::as_str));
    table.finish()
}
