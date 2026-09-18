//! CLI projection over the native bounded policy-audit engine.

use crate::AuditArguments;
use crate::commands::open;
use crate::exit::Exit;
use crate::report::{Context, Json, Table};
use std::collections::BTreeSet;

pub(crate) fn audit_selection(args: &AuditArguments, context: Context) -> Exit {
    if args.assembly_values.is_empty()
        && args.model_values.is_empty()
        && args.altloc_values.is_empty()
        && args.identifier_values.is_empty()
        && args.missing_atom_values.is_empty()
        && args.hydrogen_values.is_empty()
        && args.atom_equivalence_values.is_empty()
        && args.symmetry_values.is_empty()
        && args.alignment_values.is_empty()
        && args.precision_values.is_empty()
        && args.periodic_values.is_empty()
        && args.vdw_radii_values.is_empty()
        && args.contact_values.is_empty()
        && args.float_tolerance_values.is_empty()
    {
        eprintln!("audit requires at least one explicit policy-value dimension");
        return Exit::Usage;
    }
    let structure = match open(&args.input, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let space = match policy_space(args, context.policy) {
        Ok(space) => space,
        Err(error) => {
            eprintln!("audit policy values are invalid: {error}");
            return Exit::Policy;
        }
    };
    let plan = match space.plan() {
        Ok(plan) => plan,
        Err(error) => {
            eprintln!("audit plan refused: {error}");
            return Exit::Resource;
        }
    };
    let result = molframe::audit(
        &plan,
        |policy| {
            crate::commands::select_text(&structure, &args.query, policy, context.execution)
                .map(|evaluation| evaluation.selection)
        },
        |selection| selection.iter().collect::<BTreeSet<_>>(),
    );
    let report = match result {
        Ok(report) => report,
        Err(findings) => {
            context.findings(&findings, &args.input.display().to_string());
            return Exit::of(&findings);
        }
    };
    emit(&report, context);
    Exit::Success
}

fn policy_space(
    args: &AuditArguments,
    baseline: &molframe::AnalysisPolicy,
) -> Result<molframe::PolicySpace, molframe::PolicyConfigError> {
    let space = molframe::PolicySpace::new(baseline.clone()).with_max_runs(args.max_runs);
    let space = identity_space(args, baseline, space)?;
    let space = chemistry_space(args, baseline, space)?;
    numeric_space(args, baseline, space)
}

fn identity_space(
    args: &AuditArguments,
    baseline: &molframe::AnalysisPolicy,
    mut space: molframe::PolicySpace,
) -> Result<molframe::PolicySpace, molframe::PolicyConfigError> {
    if !args.assembly_values.is_empty() {
        let values = args
            .assembly_values
            .iter()
            .map(|value| parse_override("assembly", value, baseline).map(|policy| policy.assembly))
            .collect::<Result<Vec<_>, _>>()?;
        space = space.vary(molframe::PolicyDimension::assembly(values));
    }
    if !args.model_values.is_empty() {
        let values = args
            .model_values
            .iter()
            .map(|value| parse_override("model", value, baseline).map(|policy| policy.model))
            .collect::<Result<Vec<_>, _>>()?;
        space = space.vary(molframe::PolicyDimension::model(values));
    }
    if !args.altloc_values.is_empty() {
        let values = args
            .altloc_values
            .iter()
            .map(|value| parse_override("altloc", value, baseline).map(|policy| policy.altloc))
            .collect::<Result<Vec<_>, _>>()?;
        space = space.vary(molframe::PolicyDimension::altloc(values));
    }
    if !args.identifier_values.is_empty() {
        let values = args
            .identifier_values
            .iter()
            .map(|value| {
                parse_override("identifiers", value, baseline).map(|policy| policy.identifiers)
            })
            .collect::<Result<Vec<_>, _>>()?;
        space = space.vary(molframe::PolicyDimension::identifiers(values));
    }
    Ok(space)
}

fn chemistry_space(
    args: &AuditArguments,
    baseline: &molframe::AnalysisPolicy,
    mut space: molframe::PolicySpace,
) -> Result<molframe::PolicySpace, molframe::PolicyConfigError> {
    if !args.missing_atom_values.is_empty() {
        let values = policies("missing_atoms", &args.missing_atom_values, baseline)?
            .into_iter()
            .map(|policy| policy.missing_atoms)
            .collect::<Vec<_>>();
        space = space.vary(molframe::PolicyDimension::missing_atoms(values));
    }
    if !args.hydrogen_values.is_empty() {
        let values = policies("hydrogens", &args.hydrogen_values, baseline)?
            .into_iter()
            .map(|policy| policy.hydrogens)
            .collect::<Vec<_>>();
        space = space.vary(molframe::PolicyDimension::hydrogens(values));
    }
    if !args.atom_equivalence_values.is_empty() {
        let values = policies("atom_equivalence", &args.atom_equivalence_values, baseline)?
            .into_iter()
            .map(|policy| policy.atom_equivalence)
            .collect::<Vec<_>>();
        space = space.vary(molframe::PolicyDimension::atom_equivalence(values));
    }
    if !args.symmetry_values.is_empty() {
        let values = policies("symmetry", &args.symmetry_values, baseline)?
            .into_iter()
            .map(|policy| policy.symmetry)
            .collect::<Vec<_>>();
        space = space.vary(molframe::PolicyDimension::symmetry(values));
    }
    Ok(space)
}

fn numeric_space(
    args: &AuditArguments,
    baseline: &molframe::AnalysisPolicy,
    mut space: molframe::PolicySpace,
) -> Result<molframe::PolicySpace, molframe::PolicyConfigError> {
    if !args.alignment_values.is_empty() {
        let values = policies("alignment", &args.alignment_values, baseline)?
            .into_iter()
            .map(|policy| policy.alignment)
            .collect::<Vec<_>>();
        space = space.vary(molframe::PolicyDimension::alignment(values));
    }
    if !args.precision_values.is_empty() {
        let values = policies("precision", &args.precision_values, baseline)?
            .into_iter()
            .map(|policy| policy.precision)
            .collect::<Vec<_>>();
        space = space.vary(molframe::PolicyDimension::precision(values));
    }
    if !args.periodic_values.is_empty() {
        let values = policies("periodic", &args.periodic_values, baseline)?
            .into_iter()
            .map(|policy| policy.periodic)
            .collect::<Vec<_>>();
        space = space.vary(molframe::PolicyDimension::periodic(values));
    }
    if !args.vdw_radii_values.is_empty() {
        let values = policies("vdw_radii", &args.vdw_radii_values, baseline)?
            .into_iter()
            .map(|policy| policy.vdw_radii)
            .collect::<Vec<_>>();
        space = space.vary(molframe::PolicyDimension::vdw_radii(values));
    }
    if !args.contact_values.is_empty() {
        let values = policies("contact_def", &args.contact_values, baseline)?
            .into_iter()
            .map(|policy| policy.contact_def)
            .collect::<Vec<_>>();
        space = space.vary(molframe::PolicyDimension::contact_def(values));
    }
    if !args.float_tolerance_values.is_empty() {
        let values = args
            .float_tolerance_values
            .iter()
            .map(|value| parse_tolerance(value, baseline))
            .collect::<Result<Vec<_>, _>>()?;
        space = space.vary(molframe::PolicyDimension::float_tolerance(values));
    }
    Ok(space)
}

fn policies(
    field: &'static str,
    values: &[String],
    baseline: &molframe::AnalysisPolicy,
) -> Result<Vec<molframe::AnalysisPolicy>, molframe::PolicyConfigError> {
    values
        .iter()
        .map(|value| parse_override(field, value, baseline))
        .collect()
}

fn parse_override(
    field: &'static str,
    value: &str,
    baseline: &molframe::AnalysisPolicy,
) -> Result<molframe::AnalysisPolicy, molframe::PolicyConfigError> {
    let overrides = match field {
        "assembly" => molframe::PolicyOverrides {
            assembly: Some(value.to_owned()),
            ..molframe::PolicyOverrides::default()
        },
        "model" => molframe::PolicyOverrides {
            model: Some(value.to_owned()),
            ..molframe::PolicyOverrides::default()
        },
        "altloc" => molframe::PolicyOverrides {
            altloc: Some(value.to_owned()),
            ..molframe::PolicyOverrides::default()
        },
        "identifiers" => molframe::PolicyOverrides {
            identifiers: Some(value.to_owned()),
            ..molframe::PolicyOverrides::default()
        },
        "missing_atoms" => molframe::PolicyOverrides {
            missing_atoms: Some(value.to_owned()),
            ..molframe::PolicyOverrides::default()
        },
        "hydrogens" => molframe::PolicyOverrides {
            hydrogens: Some(value.to_owned()),
            ..molframe::PolicyOverrides::default()
        },
        "atom_equivalence" => molframe::PolicyOverrides {
            atom_equivalence: Some(value.to_owned()),
            ..molframe::PolicyOverrides::default()
        },
        "symmetry" => molframe::PolicyOverrides {
            symmetry: Some(value.to_owned()),
            ..molframe::PolicyOverrides::default()
        },
        "alignment" => molframe::PolicyOverrides {
            alignment: Some(value.to_owned()),
            ..molframe::PolicyOverrides::default()
        },
        "precision" => molframe::PolicyOverrides {
            precision: Some(value.to_owned()),
            ..molframe::PolicyOverrides::default()
        },
        "periodic" => molframe::PolicyOverrides {
            periodic: Some(value.to_owned()),
            ..molframe::PolicyOverrides::default()
        },
        "vdw_radii" => molframe::PolicyOverrides {
            vdw_radii: Some(value.to_owned()),
            ..molframe::PolicyOverrides::default()
        },
        "contact_def" => molframe::PolicyOverrides {
            contact_def: Some(value.to_owned()),
            ..molframe::PolicyOverrides::default()
        },
        _ => {
            return Err(molframe::PolicyConfigError::InvalidValue {
                field,
                value: value.to_owned(),
            });
        }
    };
    overrides.apply_to(baseline.clone())
}

fn parse_tolerance(
    value: &str,
    baseline: &molframe::AnalysisPolicy,
) -> Result<molframe::Tolerance, molframe::PolicyConfigError> {
    let Some((relative, absolute)) = value.split_once(':') else {
        return Err(molframe::PolicyConfigError::InvalidValue {
            field: "float_tolerance",
            value: value.to_owned(),
        });
    };
    let relative =
        relative
            .parse::<f64>()
            .map_err(|_| molframe::PolicyConfigError::InvalidValue {
                field: "float_tolerance",
                value: value.to_owned(),
            })?;
    let absolute =
        absolute
            .parse::<f64>()
            .map_err(|_| molframe::PolicyConfigError::InvalidValue {
                field: "float_tolerance",
                value: value.to_owned(),
            })?;
    let policy = molframe::PolicyOverrides {
        float_tolerance_relative: Some(relative),
        float_tolerance_absolute: Some(absolute),
        ..molframe::PolicyOverrides::default()
    }
    .apply_to(baseline.clone())?;
    Ok(policy.float_tolerance)
}

fn emit<R>(report: &molframe::AuditReport<R, u32>, context: Context) {
    if context.is_json() {
        let mut json = Json::new();
        let dimensions = report
            .dimensions
            .iter()
            .map(|dimension| {
                format!(
                    "{{\"field\":\"{}\",\"mean_change\":{},\"sensitive_items\":{}}}",
                    dimension.field.name(),
                    dimension.mean_change,
                    dimension.sensitive_items.len()
                )
            })
            .collect::<Vec<_>>();
        json.number("runs", report.runs.len())
            .number("stability", report.stability)
            .number("sensitive_items", report.sensitive_items.len())
            .raw("dimensions", &format!("[{}]", dimensions.join(",")));
        context.result(&json.finish());
    } else {
        let delimiter = match context.delimiter() {
            Some(delimiter) => delimiter,
            None => '\t',
        };
        let mut table = Table::new(
            delimiter,
            &["field", "mean_change", "sensitive_items", "runs"],
        );
        for dimension in &report.dimensions {
            let values = [
                dimension.field.name().to_owned(),
                dimension.mean_change.to_string(),
                dimension.sensitive_items.len().to_string(),
                report.runs.len().to_string(),
            ];
            table.row(values.iter().map(String::as_str));
        }
        context.result(&table.finish());
    }
}
