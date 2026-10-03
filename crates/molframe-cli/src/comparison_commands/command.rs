//! Thin rendering over native correspondence-based comparison kernels.

use super::render::{emit_mappings, emit_metrics, emit_rmsd};
use crate::MetricChoice;
use crate::commands::open;
use crate::exit::Exit;
use crate::report::{Context, Json};
use molframe::sequence::Scoring;
use std::path::Path;

pub(crate) fn superpose(
    mobile: &Path,
    reference: &Path,
    output: &Path,
    query: &str,
    options: molframe::geometry::SuperposeOptions,
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
    let mobile_selection = match crate::commands::select_text(
        &mobile_structure,
        query,
        context.policy,
        context.execution,
    ) {
        Ok(evaluation) => evaluation,
        Err(findings) => {
            context.findings(&findings, &mobile.display().to_string());
            return Exit::of(&findings);
        }
    };
    let reference_selection = match crate::commands::select_text(
        &reference_structure,
        query,
        context.policy,
        context.execution,
    ) {
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
    let fit = match molframe::geometry::superpose_with_options(
        &mobile_points,
        &reference_points,
        options,
    ) {
        Ok(fit) => fit,
        Err(error) => {
            eprintln!("superposition failed: {}", superpose_error(error));
            return Exit::Consistency;
        }
    };
    let all = molframe_core::selection::AtomSelection::All(mobile_structure.atom_count());
    let aligned = match molframe::transform(&mobile_structure, &all, &fit.transform) {
        Ok(aligned) => aligned,
        Err(findings) => {
            context.findings(&findings, &mobile.display().to_string());
            return Exit::of(&findings);
        }
    };
    match molframe::write(output, &aligned) {
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

const fn superpose_error(error: molframe::geometry::SuperposeError) -> &'static str {
    match error {
        molframe::geometry::SuperposeError::LengthMismatch => {
            "selections have different atom counts"
        }
        molframe::geometry::SuperposeError::TooFewPoints => "selection has fewer than three atoms",
        molframe::geometry::SuperposeError::Degenerate => "selected coordinates are collinear",
        molframe::geometry::SuperposeError::InvalidOptions => "superposition options are invalid",
        molframe::geometry::SuperposeError::TooManyPoints => {
            "selection exceeds the supported atom count"
        }
        molframe::geometry::SuperposeError::Eigen(_) => "quaternion eigensolver failed",
    }
}

fn selected_positions(
    structure: &molframe::Structure,
    selection: &molframe_core::selection::AtomSelection,
) -> Vec<[f32; 3]> {
    selection
        .iter()
        .map(|atom| structure.coordinates()[atom as usize])
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
        molframe::geometry::rmsd(&moving, &fixed).map(|value| (value, false))
    } else {
        molframe::geometry::superpose(&moving, &fixed).map(|fit| (fit.rmsd, true))
    };
    let Ok((value, fitted)) = measured else {
        eprintln!("the selected coordinates cannot be compared");
        return Exit::Consistency;
    };
    emit_rmsd(value, moving.len(), fitted, context);
    Exit::Success
}

fn comparison_points(
    mobile: &molframe::Structure,
    reference: &molframe::Structure,
    mobile_path: &Path,
    reference_path: &Path,
    query: Option<&str>,
    context: Context,
) -> Result<ComparisonPoints, Exit> {
    let Some(query) = query else {
        return Ok((
            mobile.coordinates().to_vec(),
            reference.coordinates().to_vec(),
        ));
    };
    let moving = crate::commands::select_text(mobile, query, context.policy, context.execution)
        .map_err(|findings| {
            context.findings(&findings, &mobile_path.display().to_string());
            Exit::of(&findings)
        })?;
    let fixed = crate::commands::select_text(reference, query, context.policy, context.execution)
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

type ComparisonPoints = (Vec<[f32; 3]>, Vec<[f32; 3]>);

#[derive(Clone, Copy)]
pub(crate) struct ComparisonOptions<'a> {
    pub(crate) metrics: &'a [MetricChoice],
    pub(crate) lddt_radius: Option<f32>,
    pub(crate) lddt_minimum_distance: Option<f64>,
    pub(crate) lddt_tolerances: &'a [f64],
    pub(crate) lddt_empty: Option<molframe::compare::EmptyLddtPolicy>,
    pub(crate) receptor: Option<&'a str>,
    pub(crate) ligand: Option<&'a str>,
    pub(crate) contact_distance: Option<f32>,
    pub(crate) ligand_scale: Option<f64>,
    pub(crate) interface_scale: Option<f64>,
    pub(crate) mapped: Option<MappedScoring<'a>>,
    pub(crate) extra: &'a crate::ExtraMetricArguments,
    pub(crate) chemistry: (Option<&'a Path>, Option<&'a str>),
}

/// How to put chains and atoms in correspondence before scoring.
#[derive(Clone, Copy)]
pub(crate) struct MappedScoring<'a> {
    pub(crate) ccd: &'a Path,
    pub(crate) ccd_version: &'a str,
    pub(crate) min_identity: f64,
    pub(crate) scoring: Scoring,
    pub(crate) automorphism_limit: usize,
}

pub(crate) fn compare(
    model: &Path,
    reference: &Path,
    options: ComparisonOptions<'_>,
    context: Context,
) -> Exit {
    if let Err(exit) = check_controls_have_metrics(&options) {
        return exit;
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
        let measured = match measure(
            *metric,
            &model,
            &reference,
            &options,
            context,
            context.execution,
        ) {
            Ok(value) => value,
            Err(exit) => return exit,
        };
        rows.extend(measured);
    }
    emit_metrics(context, &rows);
    Exit::Success
}

/// A control for a metric that was not selected is refused, not ignored.
fn check_controls_have_metrics(options: &ComparisonOptions<'_>) -> Result<(), Exit> {
    let has =
        |wanted: &[MetricChoice]| options.metrics.iter().any(|metric| wanted.contains(metric));
    let chain_pair = has(&[MetricChoice::DockQ, MetricChoice::Qs]);
    let cutoff = chain_pair || has(&[MetricChoice::ContactSimilarity]);
    let refused = [
        (
            (options.ligand_scale.is_some() || options.interface_scale.is_some())
                && !has(&[MetricChoice::DockQ]),
            "DockQ scales require `--metrics dockq`",
        ),
        (
            (options.receptor.is_some() || options.ligand.is_some()) && !chain_pair,
            "--receptor and --ligand require `--metrics dockq` or `qs`",
        ),
        (
            options.contact_distance.is_some() && !cutoff,
            "--contact-distance requires `dockq`, `qs` or `contact-similarity`",
        ),
    ];
    for (violated, message) in refused {
        if violated {
            eprintln!("{message}");
            return Err(Exit::Usage);
        }
    }
    Ok(())
}

fn measure(
    metric: MetricChoice,
    model: &molframe::Structure,
    reference: &molframe::Structure,
    options: &ComparisonOptions<'_>,
    context: Context,
    execution: &molframe_core::ExecutionContext,
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
            molframe::compare::lddt_with_options(
                model.coordinates(),
                reference.coordinates(),
                &molframe::compare::LddtOptions {
                    inclusion_radius: f64::from(radius),
                    minimum_reference_distance: minimum,
                    tolerances: options.lddt_tolerances.into(),
                    empty_policy: empty,
                },
                execution,
            )
        }
        MetricChoice::TmScore => {
            molframe::compare::tm_score(model.coordinates(), reference.coordinates())
        }
        MetricChoice::GdtTs => {
            molframe::compare::gdt_ts(model.coordinates(), reference.coordinates())
        }
        MetricChoice::GdtHa => {
            molframe::compare::gdt_ha(model.coordinates(), reference.coordinates())
        }
        MetricChoice::DockQ => {
            return super::docking::measure_dockq(model, reference, options, context);
        }
        MetricChoice::Qs => return super::overlap::qs(model, reference, options, context),
        MetricChoice::ContactSimilarity => {
            return super::overlap::contact_similarity(model, reference, options, context);
        }
        MetricChoice::Cad => return super::overlap::cad(model, reference, options, context),
        MetricChoice::Ce => return super::ce::measure(model, reference, options, context),
        MetricChoice::InterfaceRmsd | MetricChoice::PocketRmsd => {
            return super::regions::region_rmsd(metric, model, reference, options, context);
        }
        MetricChoice::LigandRmsd => {
            return super::regions::ligand_rmsd(model, reference, options, context);
        }
    };
    result
        .map(|value| vec![(metric_name(metric), value)])
        .map_err(|error| {
            eprintln!("comparison failed: {error}");
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
    let assignment = match molframe::compare::assign_chains(
        reference.engine(),
        model.engine(),
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

pub(super) const fn metric_name(metric: MetricChoice) -> &'static str {
    match metric {
        MetricChoice::Lddt => "lddt",
        MetricChoice::TmScore => "tm_score",
        MetricChoice::GdtTs => "gdt_ts",
        MetricChoice::GdtHa => "gdt_ha",
        MetricChoice::DockQ => "dockq",
        MetricChoice::Qs => "qs",
        MetricChoice::Cad => "cad",
        MetricChoice::Ce => "ce_rmsd",
        MetricChoice::ContactSimilarity => "contact_similarity",
        MetricChoice::LigandRmsd => "ligand_rmsd",
        MetricChoice::InterfaceRmsd => "interface_rmsd",
        MetricChoice::PocketRmsd => "pocket_rmsd",
    }
}
