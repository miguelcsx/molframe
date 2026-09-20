//! CLI projection over native declarative functional-geometry evaluation.

use crate::FxCommand;
use crate::commands::open;
use crate::exit::Exit;
use crate::report::{Context, Json, Table};

pub(crate) fn execute(command: FxCommand, context: Context) -> Exit {
    match command {
        FxCommand::Evaluate {
            input,
            motif,
            chemistry,
            mapping_limit,
            measurement_limit,
            plane_relative_tolerance,
            plane_maximum_sweeps,
        } => evaluate(
            &input,
            &motif,
            chemistry.ccd.as_deref(),
            chemistry.ccd_version.as_deref(),
            EvaluationOptions {
                mapping_limit,
                measurement: molframe::motif::MeasurementOptions {
                    maximum_alternatives: measurement_limit,
                    plane_fit: molframe::geometry::EigenOptions {
                        relative_tolerance: plane_relative_tolerance,
                        maximum_sweeps: plane_maximum_sweeps,
                    },
                },
            },
            context,
        ),
    }
}

#[derive(Clone, Copy)]
struct EvaluationOptions {
    mapping_limit: usize,
    measurement: molframe::motif::MeasurementOptions,
}

fn evaluate(
    input: &std::path::Path,
    motif_path: &std::path::Path,
    ccd_path: Option<&std::path::Path>,
    ccd_version: Option<&str>,
    options: EvaluationOptions,
    context: Context,
) -> Exit {
    if options.mapping_limit == 0 || options.measurement.maximum_alternatives == 0 {
        eprintln!("functional evaluation limits must be positive");
        return Exit::Usage;
    }
    let structure = match open(input, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let source = match crate::chemistry::resolve_ccd(ccd_path, ccd_version, context) {
        Ok(source) => source,
        Err(exit) => return exit,
    };
    let context = source.context;
    let provider = match crate::chemistry::load_ccd(source.path, source.version, context) {
        Ok(provider) => provider,
        Err(exit) => return exit,
    };
    let specification = match molframe::motif::read_evaluation_specification(motif_path) {
        Ok(specification) => specification,
        Err(error) => {
            eprintln!("functional specification failed: {error}");
            return Exit::Policy;
        }
    };
    let report = match molframe::motif::evaluate_motif(
        structure.engine(),
        &specification.motif,
        Some(&provider),
        context.policy,
        &specification.profile,
        options.mapping_limit,
        options.measurement,
    ) {
        Ok(report) => report,
        Err(error) => {
            eprintln!("functional evaluation failed: {error}");
            return Exit::Indeterminate;
        }
    };
    emit(&report, context);
    if report.evaluations.is_empty()
        || report.evaluations.iter().any(|evaluation| {
            evaluation.verdict.status == molframe::motif::VerdictStatus::Indeterminate
        })
    {
        Exit::Indeterminate
    } else {
        Exit::Success
    }
}

fn emit(report: &molframe::motif::EvaluationReport, context: Context) {
    let rows = report
        .evaluations
        .iter()
        .flat_map(|evaluation| {
            evaluation
                .measurements
                .constraints
                .iter()
                .map(move |measurement| {
                    vec![
                        evaluation.mapping_index.to_string(),
                        measurement.name.to_string(),
                        measurement
                            .value
                            .map_or_else(|| "null".to_owned(), |value| value.numeric().to_string()),
                        measurement
                            .deviation
                            .map_or_else(|| "null".to_owned(), |value| value.to_string()),
                        measurement
                            .satisfied
                            .map_or_else(|| "null".to_owned(), |value| value.to_string()),
                        verdict_name(evaluation.verdict.status).to_owned(),
                    ]
                })
        })
        .collect::<Vec<_>>();
    let header = [
        "mapping",
        "constraint",
        "value",
        "deviation",
        "satisfied",
        "verdict",
    ];
    if context.is_json() {
        let objects = rows
            .iter()
            .map(|row| {
                let mut json = Json::new();
                for (field, value) in header.iter().zip(row) {
                    json.text(field, value);
                }
                json.finish()
            })
            .collect::<Vec<_>>();
        let mut envelope = Json::new();
        envelope
            .raw(
                "mapping_ambiguous",
                if report.mapping_ambiguous {
                    "true"
                } else {
                    "false"
                },
            )
            .raw("evaluations", &format!("[{}]", objects.join(",")));
        context.result(&envelope.finish());
    } else {
        let delimiter = match context.delimiter() {
            Some(delimiter) => delimiter,
            None => '\t',
        };
        let mut table = Table::new(delimiter, &header);
        for row in &rows {
            table.row(row.iter().map(String::as_str));
        }
        context.result(&table.finish());
    }
}

const fn verdict_name(status: molframe::motif::VerdictStatus) -> &'static str {
    match status {
        molframe::motif::VerdictStatus::Pass => "pass",
        molframe::motif::VerdictStatus::Fail => "fail",
        molframe::motif::VerdictStatus::Indeterminate => "indeterminate",
    }
}
