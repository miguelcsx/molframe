//! Thin rendering over native correspondence-based comparison kernels.

use crate::MetricChoice;
use crate::commands::open;
use crate::exit::Exit;
use crate::report::{Context, Json, Table};
use pdbiox::QueryStructure as _;
use pdbiox::seq::Scoring;
use std::path::Path;

pub(crate) fn superpose(
    mobile: &Path,
    reference: &Path,
    output: &Path,
    query: &str,
    options: pdbiox::SuperposeOptions,
    context: Context,
) -> Exit {
    use pdbiox::QueryStructure as _;
    let mobile_structure = match open(mobile, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let reference_structure = match open(reference, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let mobile_selection =
        match mobile_structure.select_text(query, context.policy, context.execution) {
            Ok(evaluation) => evaluation,
            Err(findings) => {
                context.findings(&findings, &mobile.display().to_string());
                return Exit::of(&findings);
            }
        };
    let reference_selection =
        match reference_structure.select_text(query, context.policy, context.execution) {
            Ok(evaluation) => evaluation,
            Err(findings) => {
                context.findings(&findings, &reference.display().to_string());
                return Exit::of(&findings);
            }
        };
    context.findings(&mobile_selection.warnings, &mobile.display().to_string());
    context.findings(
        &reference_selection.warnings,
        &reference.display().to_string(),
    );
    let mobile_points = selected_positions(&mobile_structure, &mobile_selection.selection);
    let reference_points = selected_positions(&reference_structure, &reference_selection.selection);
    let fit = match pdbiox::superpose_with_options(&mobile_points, &reference_points, options) {
        Ok(fit) => fit,
        Err(error) => {
            eprintln!("superposition failed: {}", superpose_error(error));
            return Exit::Consistency;
        }
    };
    let all = pdbiox::AtomSelection::All(mobile_structure.atom_count());
    let aligned = match pdbiox::transform(&mobile_structure, &all, &fit.transform) {
        Ok(aligned) => aligned,
        Err(findings) => {
            context.findings(&findings, &mobile.display().to_string());
            return Exit::of(&findings);
        }
    };
    match pdbiox::write(output, &aligned) {
        Ok(()) => {
            if context.is_json() {
                let mut json = Json::new();
                json.text("output", &output.display().to_string())
                    .number("rmsd_angstrom", fit.rmsd)
                    .number("matched_atoms", mobile_points.len());
                context.result(&json.finish());
            } else {
                context.result(&format!(
                    "wrote {} after fitting {} atoms (RMSD {} angstrom)",
                    output.display(),
                    mobile_points.len(),
                    fit.rmsd
                ));
            }
            Exit::Success
        }
        Err(findings) => {
            context.findings(&findings, &output.display().to_string());
            Exit::of(&findings)
        }
    }
}

const fn superpose_error(error: pdbiox::SuperposeError) -> &'static str {
    match error {
        pdbiox::SuperposeError::LengthMismatch => "selections have different atom counts",
        pdbiox::SuperposeError::TooFewPoints => "selection has fewer than three atoms",
        pdbiox::SuperposeError::Degenerate => "selected coordinates are collinear",
        pdbiox::SuperposeError::InvalidOptions => "superposition options are invalid",
        pdbiox::SuperposeError::TooManyPoints => "selection exceeds the supported atom count",
        pdbiox::SuperposeError::Eigen(_) => "quaternion eigensolver failed",
    }
}

fn selected_positions(
    structure: &pdbiox::Structure,
    selection: &pdbiox::AtomSelection,
) -> Vec<[f32; 3]> {
    selection
        .iter()
        .map(|atom| structure.positions()[atom as usize])
        .collect()
}

pub(crate) fn rmsd(
    mobile: &Path,
    reference: &Path,
    no_fit: bool,
    query: Option<&str>,
    context: Context,
) -> Exit {
    let mobile_structure = match open(mobile, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let reference_structure = match open(reference, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let (moving, fixed) = match comparison_points(
        &mobile_structure,
        &reference_structure,
        mobile,
        reference,
        query,
        context,
    ) {
        Ok(points) => points,
        Err(exit) => return exit,
    };
    if moving.len() != fixed.len() {
        eprintln!(
            "the selected structures hold {} and {} atoms; they must correspond one to one",
            moving.len(),
            fixed.len()
        );
        return Exit::Consistency;
    }
    let measured = if no_fit {
        pdbiox::rmsd(&moving, &fixed).map(|value| (value, false))
    } else {
        pdbiox::superpose(&moving, &fixed).map(|fit| (fit.rmsd, true))
    };
    let Ok((value, fitted)) = measured else {
        eprintln!("the selected coordinates cannot be compared");
        return Exit::Consistency;
    };
    emit_rmsd(value, moving.len(), fitted, context);
    Exit::Success
}

fn comparison_points(
    mobile: &pdbiox::Structure,
    reference: &pdbiox::Structure,
    mobile_path: &Path,
    reference_path: &Path,
    query: Option<&str>,
    context: Context,
) -> Result<ComparisonPoints, Exit> {
    let Some(query) = query else {
        return Ok((mobile.positions().to_vec(), reference.positions().to_vec()));
    };
    let moving = mobile
        .select_text(query, context.policy, context.execution)
        .map_err(|findings| {
            context.findings(&findings, &mobile_path.display().to_string());
            Exit::of(&findings)
        })?;
    let fixed = reference
        .select_text(query, context.policy, context.execution)
        .map_err(|findings| {
            context.findings(&findings, &reference_path.display().to_string());
            Exit::of(&findings)
        })?;
    context.findings(&moving.warnings, &mobile_path.display().to_string());
    context.findings(&fixed.warnings, &reference_path.display().to_string());
    Ok((
        selected_positions(mobile, &moving.selection),
        selected_positions(reference, &fixed.selection),
    ))
}

fn emit_rmsd(value: f64, atoms: usize, fitted: bool, context: Context) {
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

type ComparisonPoints = (Vec<[f32; 3]>, Vec<[f32; 3]>);

#[derive(Clone, Copy)]
pub(crate) struct ComparisonOptions<'a> {
    pub(crate) metrics: &'a [MetricChoice],
    pub(crate) lddt_radius: Option<f32>,
    pub(crate) lddt_minimum_distance: Option<f64>,
    pub(crate) lddt_tolerances: &'a [f64],
    pub(crate) lddt_empty: Option<pdbiox::compare::EmptyLddtPolicy>,
    pub(crate) receptor: Option<&'a str>,
    pub(crate) ligand: Option<&'a str>,
    pub(crate) contact_distance: Option<f32>,
    pub(crate) ligand_scale: Option<f64>,
    pub(crate) interface_scale: Option<f64>,
}

pub(crate) fn compare(
    model: &Path,
    reference: &Path,
    options: ComparisonOptions<'_>,
    context: Context,
) -> Exit {
    let includes_dockq = options
        .metrics
        .iter()
        .any(|metric| matches!(metric, MetricChoice::DockQ));
    let has_dockq_controls = options.receptor.is_some()
        || options.ligand.is_some()
        || options.contact_distance.is_some()
        || options.ligand_scale.is_some()
        || options.interface_scale.is_some();
    if has_dockq_controls && !includes_dockq {
        eprintln!("DockQ controls require `--metrics dockq`");
        return Exit::Usage;
    }
    let model = match open(model, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let reference = match open(reference, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let mut rows = Vec::new();
    for metric in options.metrics {
        let measured = match measure(*metric, &model, &reference, &options, context.execution) {
            Ok(value) => value,
            Err(exit) => return exit,
        };
        rows.extend(measured);
    }
    emit_metrics(context, &rows);
    Exit::Success
}

fn measure(
    metric: MetricChoice,
    model: &pdbiox::Structure,
    reference: &pdbiox::Structure,
    options: &ComparisonOptions<'_>,
    execution: &pdbiox::core::ExecutionContext,
) -> Result<Vec<(&'static str, f64)>, Exit> {
    let result = match metric {
        MetricChoice::Lddt => {
            let (Some(radius), Some(minimum), Some(empty)) = (
                options.lddt_radius,
                options.lddt_minimum_distance,
                options.lddt_empty,
            ) else {
                eprintln!("lDDT requires --lddt-radius, --lddt-minimum-distance and --lddt-empty");
                return Err(Exit::Usage);
            };
            if options.lddt_tolerances.is_empty() {
                eprintln!("lDDT requires one or more --lddt-tolerances");
                return Err(Exit::Usage);
            }
            pdbiox::compare::lddt_with_options(
                model.positions(),
                reference.positions(),
                &pdbiox::compare::LddtOptions {
                    inclusion_radius: f64::from(radius),
                    minimum_reference_distance: minimum,
                    tolerances: options.lddt_tolerances.into(),
                    empty_policy: empty,
                },
                execution,
            )
        }
        MetricChoice::TmScore => {
            pdbiox::compare::tm_score(model.positions(), reference.positions())
        }
        MetricChoice::GdtTs => pdbiox::compare::gdt_ts(model.positions(), reference.positions()),
        MetricChoice::GdtHa => pdbiox::compare::gdt_ha(model.positions(), reference.positions()),
        MetricChoice::DockQ => return measure_dockq(model, reference, options),
    };
    result
        .map(|value| vec![(metric_name(metric), value)])
        .map_err(|error| {
            eprintln!("comparison failed: {error}");
            Exit::Consistency
        })
}

fn measure_dockq(
    model: &pdbiox::Structure,
    reference: &pdbiox::Structure,
    options: &ComparisonOptions<'_>,
) -> Result<Vec<(&'static str, f64)>, Exit> {
    let Some(receptor) = options.receptor else {
        eprintln!("--receptor is required when dock-q is selected");
        return Err(Exit::Usage);
    };
    let Some(ligand) = options.ligand else {
        eprintln!("--ligand is required when dock-q is selected");
        return Err(Exit::Usage);
    };
    let (Some(contact_distance), Some(ligand_scale), Some(interface_scale)) = (
        options.contact_distance,
        options.ligand_scale,
        options.interface_scale,
    ) else {
        eprintln!(
            "--contact-distance, --ligand-scale and --interface-scale are required when dock-q is selected"
        );
        return Err(Exit::Usage);
    };
    pdbiox::compare::dockq(
        model,
        reference,
        receptor,
        ligand,
        pdbiox::compare::DockQOptions {
            contact_distance,
            ligand_scale,
            interface_scale,
        },
    )
    .map(|result| {
        vec![
            ("dockq", result.score),
            ("dockq_fnat", result.fnat),
            ("dockq_ligand_rmsd", result.ligand_rmsd),
            ("dockq_interface_rmsd", result.interface_rmsd),
        ]
    })
    .map_err(|error| {
        eprintln!("DockQ failed: {error}");
        Exit::Consistency
    })
}

pub(crate) fn map_chains(
    model: &Path,
    reference: &Path,
    ccd: &Path,
    ccd_version: &str,
    min_identity: f64,
    scoring: Scoring,
    context: Context,
) -> Exit {
    if !min_identity.is_finite() || !(0.0..=1.0).contains(&min_identity) {
        eprintln!("--min-identity must be between zero and one");
        return Exit::Usage;
    }
    let reference_origin = reference.display().to_string();
    let model = match open(model, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let reference = match open(reference, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    let provider = match crate::chemistry::load_ccd(ccd, ccd_version, context) {
        Ok(provider) => provider,
        Err(exit) => return exit,
    };
    let assignment = match pdbiox::compare::assign_chains(
        &reference,
        &model,
        &provider,
        context.policy.identifiers,
        scoring,
        min_identity,
    ) {
        Ok(assignment) => assignment,
        Err(finding) => {
            context.findings(&[finding], &reference_origin);
            return Exit::Consistency;
        }
    };
    let mut rows = assignment
        .primary
        .into_iter()
        .map(|mapping| {
            (
                "primary",
                mapping.reference,
                mapping.target,
                mapping.identity,
            )
        })
        .collect::<Vec<_>>();
    rows.extend(assignment.alternatives.into_iter().map(|mapping| {
        (
            "alternative",
            mapping.reference,
            mapping.target,
            mapping.identity,
        )
    }));
    emit_mappings(context, &rows);
    Exit::Success
}

fn emit_metrics(context: Context, rows: &[(&str, f64)]) {
    if context.is_json() {
        let mut json = Json::new();
        for (name, value) in rows {
            json.number(name, value);
        }
        context.result(&json.finish());
    } else {
        let delimiter = match context.delimiter() {
            Some(value) => value,
            None => '\t',
        };
        let mut table = Table::new(delimiter, &["metric", "value"]);
        for (name, value) in rows {
            let value = value.to_string();
            table.row([*name, &value]);
        }
        context.result(&table.finish());
    }
}

fn emit_mappings(context: Context, rows: &[(&str, String, String, f64)]) {
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
        let delimiter = match context.delimiter() {
            Some(value) => value,
            None => '\t',
        };
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

pub(super) const fn metric_name(metric: MetricChoice) -> &'static str {
    match metric {
        MetricChoice::Lddt => "lddt",
        MetricChoice::TmScore => "tm_score",
        MetricChoice::GdtTs => "gdt_ts",
        MetricChoice::GdtHa => "gdt_ha",
        MetricChoice::DockQ => "dockq",
    }
}
