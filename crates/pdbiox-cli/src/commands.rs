//! What each subcommand does.
//!
//! Every one of these reads through the library, formats what came back, and
//! returns an exit code. None of them computes anything: a convenience worth
//! having belongs in the library, where the Rust callers get it too.

use crate::exit::Exit;
use crate::report::{Context, Json, Table, json_array};
use pdbiox::{
    AnalysisPolicy, ChainRef, Format, PdbOptions, ReadOptions, Structure, validate as check,
};
use std::fmt::Write as _;
use std::path::Path;

/// Reads a file, printing whatever was wrong with it.
fn open(path: &Path, context: Context) -> Result<Structure, Exit> {
    let options = ReadOptions::new().mode(context.mode);
    let origin = path.display().to_string();
    match pdbiox::read_with_options(path, &options) {
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
    let data = structure.data();

    if context.is_json() {
        let mut object = Json::new();
        if let Some(id) = &data.entry.id {
            object.text("id", id);
        }
        object
            .number("models", structure.model_count())
            .number("chains", structure.chain_count())
            .number("entities", structure.entity_count())
            .number("residues", structure.residue_count())
            .number("atoms", structure.atom_count());
        if detail {
            object.raw("chain_detail", &chain_detail_json(&structure));
        }
        context.result(&object.finish());
        return Exit::Success;
    }
    if let Some(delimiter) = context.delimiter() {
        context.result(&info_table(&structure, detail, delimiter));
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
    let _ = writeln!(text, "entities  {}", structure.entity_count());
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
                label_of(&structure, chain),
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

fn info_table(structure: &Structure, detail: bool, delimiter: char) -> String {
    let id = match &structure.data().entry.id {
        Some(id) => id.as_ref(),
        None => "",
    };
    if detail {
        let mut table = Table::new(delimiter, &["entry", "chain", "residues", "atoms"]);
        for chain in structure.data().chains() {
            let values = [
                id.to_owned(),
                label_of(structure, chain),
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
        structure.entity_count().to_string(),
        structure.residue_count().to_string(),
        structure.atom_count().to_string(),
    ];
    table.row(values.iter().map(String::as_str));
    table.finish()
}

fn chain_detail_json(structure: &Structure) -> String {
    let entries: Vec<String> = structure
        .data()
        .chains()
        .map(|chain| {
            let mut object = Json::new();
            object
                .text("chain", &label_of(structure, chain))
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

fn label_of(structure: &Structure, chain: ChainRef<'_>) -> String {
    let symbol = chain.auth_asym_id().or_else(|| chain.label_asym_id());
    match symbol.and_then(|symbol| structure.resolve(symbol)) {
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
    context: Context,
) -> Exit {
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

    let bytes = match target {
        Format::Mmcif => pdbiox::write_mmcif(&structure).into_bytes(),
        Format::BinaryCif => match pdbiox::write_bcif(&structure) {
            Ok(bytes) => bytes,
            Err(refusals) => {
                context.findings(&refusals, &input.display().to_string());
                return Exit::of(&refusals);
            }
        },
        Format::Pdb => {
            let options = match pdb_options(chain_map, hybrid36) {
                Ok(options) => options,
                Err(exit) => return exit,
            };
            match pdbiox::write_pdb(&structure, &options) {
                Ok(text) => text.into_bytes(),
                Err(refusals) => {
                    context.findings(&refusals, &input.display().to_string());
                    return Exit::of(&refusals);
                }
            }
        }
        Format::Pqr | Format::Pdbqt => {
            let options = match pdb_options(chain_map, hybrid36) {
                Ok(options) => options,
                Err(exit) => return exit,
            };
            let rendered = match target {
                Format::Pqr => pdbiox::write_pqr(&structure, &options),
                Format::Pdbqt => pdbiox::write_pdbqt(&structure, &options),
                _ => return Exit::Usage,
            };
            match rendered {
                Ok(text) => text.into_bytes(),
                Err(refusals) => {
                    context.findings(&refusals, &input.display().to_string());
                    return Exit::of(&refusals);
                }
            }
        }
        _ => return Exit::Usage,
    };
    match std::fs::write(output, bytes) {
        Ok(()) => Exit::Success,
        Err(error) => {
            eprintln!("could not write {}: {error}", output.display());
            Exit::Failure
        }
    }
}

fn pdb_options(chain_map: &[String], hybrid36: bool) -> Result<PdbOptions, Exit> {
    let mut options = PdbOptions::new().hybrid36(hybrid36);
    for mapping in chain_map {
        let Some((from, to)) = mapping.split_once('=') else {
            eprintln!("a chain mapping must be written FROM=TO, not {mapping:?}");
            return Err(Exit::Usage);
        };
        options = options.chain_map(from, to);
    }
    Ok(options)
}

/// Checks a structure against the invariants it must satisfy.
pub fn validate(path: &Path, context: Context) -> Exit {
    let structure = match open(path, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let violations = check(structure.data());

    if context.is_json() {
        let mut object = Json::new();
        object.number("violations", violations.len()).text(
            "status",
            if violations.is_empty() {
                "sound"
            } else {
                "inconsistent"
            },
        );
        context.result(&object.finish());
    } else if let Some(delimiter) = context.delimiter() {
        let values = [
            violations.len().to_string(),
            if violations.is_empty() {
                "sound".to_owned()
            } else {
                "inconsistent".to_owned()
            },
        ];
        context.result(&one_row(delimiter, &["violations", "status"], &values));
    } else if violations.is_empty() {
        context.result("sound: every structural invariant holds");
    } else {
        context.result(&format!("inconsistent: {} violations", violations.len()));
    }

    context.findings(&violations, &path.display().to_string());
    if violations.is_empty() {
        Exit::Success
    } else {
        Exit::Consistency
    }
}

/// Measures a structure's overall shape.
pub fn measure(path: &Path, context: Context) -> Exit {
    let structure = match open(path, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let positions = structure.positions();

    let Some(centre) = pdbiox::centroid(positions) else {
        eprintln!("the structure holds no positions to measure");
        return Exit::Consistency;
    };
    let radius = pdbiox::radius_of_gyration(positions, &[]);

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

/// Compares two structures, fitting one onto the other unless told not to.
pub fn rmsd(mobile: &Path, reference: &Path, no_fit: bool, context: Context) -> Exit {
    let (Ok(mobile), Ok(reference)) = (open(mobile, context), open(reference, context)) else {
        return Exit::Input;
    };
    let (moving, fixed) = (mobile.positions(), reference.positions());

    // Fitting needs the two sets to correspond atom by atom. Deciding which
    // atom matches which is a mapping question, and answering it by position is
    // only right when the two files describe the same atoms in the same order.
    if moving.len() != fixed.len() {
        eprintln!(
            "the structures hold {} and {} atoms; they must correspond one to one",
            moving.len(),
            fixed.len()
        );
        return Exit::Consistency;
    }

    let measured = if no_fit {
        pdbiox::rmsd(moving, fixed).map(|value| (value, false))
    } else {
        pdbiox::superpose(moving, fixed).map(|fit| (fit.rmsd, true))
    };
    let Ok((value, fitted)) = measured else {
        eprintln!("the structures cannot be compared: too few points, or they disagree in size");
        return Exit::Consistency;
    };

    if context.is_json() {
        let mut object = Json::new();
        object.number("rmsd", format!("{value:.4}"));
        object.number("atoms", moving.len());
        object.text("fitted", if fitted { "yes" } else { "no" });
        context.result(&object.finish());
    } else if let Some(delimiter) = context.delimiter() {
        let values = [
            format!("{value:.4}"),
            moving.len().to_string(),
            fitted.to_string(),
        ];
        context.result(&one_row(delimiter, &["rmsd", "atoms", "fitted"], &values));
    } else {
        let how = if fitted {
            "after fitting"
        } else {
            "as they sit"
        };
        context.result(&format!(
            "rmsd {value:.3} over {} atoms, {how}",
            moving.len()
        ));
    }
    Exit::Success
}

/// Prints the policy an analysis runs under by default.
pub fn policy(context: Context) -> Exit {
    let policy = AnalysisPolicy::default();
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
