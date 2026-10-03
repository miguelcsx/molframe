//! Detailed projections over native validation kernels.
//!
//! `validate --checks a,b,c` runs each named check against one structure and
//! prints one row per finding. A check that embodies a scientific definition
//! refuses to run without its explicit parameters.

mod basic;
mod reference;
mod structural;

use crate::ValidationChoice;
use crate::args::ValidateArguments;
use crate::commands::open;
use crate::exit::Exit;
use crate::report::{Context, Json, Table};

/// What `validate` was asked to check, and with which explicit parameters.
pub(crate) type ValidationOptions<'a> = &'a ValidateArguments;

/// Whether a check reads CCD-derived roles, charges, bonds or stereochemistry.
const fn needs_dictionary(check: ValidationChoice) -> bool {
    matches!(
        check,
        ValidationChoice::Geometry
            | ValidationChoice::Chirality
            | ValidationChoice::Nucleic
            | ValidationChoice::Rotamer
            | ValidationChoice::CisPeptide
            | ValidationChoice::Ramachandran
    )
}

/// Whether the structure is annotated from the dictionary before checking.
///
/// Valence and ligand checks need only a bond table, so they annotate only
/// when the caller supplied a dictionary to take one from.
fn wants_annotation(options: ValidationOptions<'_>) -> bool {
    options.checks.iter().copied().any(needs_dictionary)
        || (options.chemistry.ccd.is_some()
            && options
                .checks
                .iter()
                .any(|check| matches!(check, ValidationChoice::Valence | ValidationChoice::Ligand)))
}

/// The inputs every check can draw on.
pub(super) struct Subject<'a> {
    pub(super) structure: &'a molframe::Structure,
    pub(super) options: ValidationOptions<'a>,
    pub(super) context: Context,
}

pub(crate) fn validate(options: ValidationOptions<'_>, context: Context) -> Exit {
    let input = options.input.as_path();
    let mut context = context;
    let mut structure = match open(input, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    if wants_annotation(options) {
        let source = match crate::chemistry::resolve_ccd(
            options.chemistry.ccd.as_deref(),
            options.chemistry.ccd_version.as_deref(),
            context,
        ) {
            Ok(source) => source,
            Err(exit) => return exit,
        };
        context = source.context;
        structure =
            match crate::chemistry::annotate(&structure, source.path, source.version, context) {
                Ok(structure) => structure,
                Err(exit) => return exit,
            };
    }
    let subject = Subject {
        structure: &structure,
        options,
        context,
    };
    let mut rows = Vec::new();
    for check in &options.checks {
        match run_check(*check, &subject) {
            Ok(mut findings) => rows.append(&mut findings),
            Err(exit) => return exit,
        }
    }
    emit(&rows, context);
    if rows.is_empty() {
        Exit::Success
    } else {
        Exit::Consistency
    }
}

fn run_check(check: ValidationChoice, subject: &Subject<'_>) -> Result<Vec<Row>, Exit> {
    let Subject {
        structure,
        options,
        context,
    } = *subject;
    match check {
        ValidationChoice::Core => Ok(basic::core_rows(
            structure,
            context,
            options.input.as_path(),
        )),
        ValidationChoice::Quality => Ok(basic::quality_rows(structure)),
        ValidationChoice::Geometry => basic::geometry_rows(structure, options),
        ValidationChoice::Clashes => basic::clash_rows(structure, options, context.execution),
        ValidationChoice::Completeness => basic::completeness_rows(structure, context),
        ValidationChoice::Altloc => basic::altloc_rows(structure, options, context),
        ValidationChoice::CcdCompleteness => basic::ccd_rows(structure, options, context),
        ValidationChoice::Bfactor => basic::b_factor_rows(structure, options),
        ValidationChoice::Valence => structural::valence_rows(subject),
        ValidationChoice::Ligand => structural::ligand_rows(subject),
        ValidationChoice::Chirality => structural::chirality_rows(subject),
        ValidationChoice::CisPeptide => structural::cis_peptide_rows(subject),
        ValidationChoice::Nucleic => structural::nucleic_rows(subject),
        ValidationChoice::Ramachandran => reference::ramachandran_rows(subject),
        ValidationChoice::Rotamer => reference::rotamer_rows(subject),
        ValidationChoice::PlaneRestraints => reference::plane_restraint_rows(subject),
        ValidationChoice::Tls => reference::tls_rows(subject),
    }
}

#[derive(Debug)]
pub(super) struct Row {
    check: String,
    item: String,
    value: String,
}

impl Row {
    pub(super) fn new(
        check: impl Into<String>,
        item: impl Into<String>,
        value: impl Into<String>,
    ) -> Self {
        Self {
            check: check.into(),
            item: item.into(),
            value: value.into(),
        }
    }
}

fn emit(rows: &[Row], context: Context) {
    if context.is_json() {
        let records: Vec<String> = rows
            .iter()
            .map(|row| {
                let mut object = Json::new();
                object
                    .text("check", &row.check)
                    .text("item", &row.item)
                    .text("value", &row.value);
                object.finish()
            })
            .collect();
        context.result(&context.json_records(&records));
    } else {
        let delimiter = context.table_delimiter();
        let mut table = Table::new(delimiter, &["check", "item", "value"]);
        for row in rows {
            table.row([row.check.as_str(), row.item.as_str(), row.value.as_str()]);
        }
        context.result(&table.finish());
    }
}
