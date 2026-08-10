//! Detailed projections over native validation kernels.

use crate::ValidationChoice;
use crate::commands::open;
use crate::exit::Exit;
use crate::report::{Context, Json, Table};
use pdbiox::{RadiusSet, SpatialBackend};
use std::path::Path;

#[derive(Clone, Copy)]
pub(crate) struct ValidationOptions<'a> {
    pub checks: &'a [ValidationChoice],
    pub ccd: Option<&'a Path>,
    pub ccd_version: Option<&'a str>,
    pub bond_tolerance: Option<f32>,
    pub planarity_tolerance: Option<f64>,
    pub plane_relative_tolerance: Option<f64>,
    pub plane_maximum_sweeps: Option<usize>,
    pub clash_tolerance: Option<f32>,
    pub radii: Option<RadiusSet>,
    pub altloc_expected_sum: Option<f64>,
    pub altloc_tolerance: Option<f64>,
    pub b_factor_z_score: Option<f64>,
}

pub(crate) fn validate(input: &Path, options: ValidationOptions<'_>, context: Context) -> Exit {
    validate_with(input, options, context)
}

fn validate_with(input: &Path, options: ValidationOptions<'_>, context: Context) -> Exit {
    let mut context = context;
    let mut structure = match open(input, context) {
        Ok(structure) => structure,
        Err(exit) => return exit,
    };
    if options.checks.contains(&ValidationChoice::Geometry) {
        let source = match crate::chemistry::resolve_ccd(options.ccd, options.ccd_version, context)
        {
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
    let mut rows = Vec::new();
    for check in options.checks {
        let mut findings = match check {
            ValidationChoice::Core => core_rows(&structure, context, input),
            ValidationChoice::Quality => quality_rows(&structure),
            ValidationChoice::Geometry => match geometry_rows(&structure, &options) {
                Ok(rows) => rows,
                Err(exit) => return exit,
            },
            ValidationChoice::Clashes => match clash_rows(&structure, &options) {
                Ok(rows) => rows,
                Err(exit) => return exit,
            },
            ValidationChoice::Completeness => match completeness_rows(&structure, context) {
                Ok(rows) => rows,
                Err(exit) => return exit,
            },
            ValidationChoice::Altloc => match altloc_rows(&structure, &options, context) {
                Ok(rows) => rows,
                Err(exit) => return exit,
            },
            ValidationChoice::CcdCompleteness => match ccd_rows(&structure, &options, context) {
                Ok(rows) => rows,
                Err(exit) => return exit,
            },
            ValidationChoice::Bfactor => match b_factor_rows(&structure, &options) {
                Ok(rows) => rows,
                Err(exit) => return exit,
            },
        };
        rows.append(&mut findings);
    }
    emit(&rows, context);
    if rows.is_empty() {
        Exit::Success
    } else {
        Exit::Consistency
    }
}

fn b_factor_rows(
    structure: &pdbiox::Structure,
    options: &ValidationOptions<'_>,
) -> Result<Vec<Row>, Exit> {
    let Some(threshold) = options.b_factor_z_score else {
        eprintln!("B-factor validation requires --b-factor-z-score");
        return Err(Exit::Usage);
    };
    let selection = pdbiox::AtomSelection::from_sorted((0..structure.atom_count()).collect());
    let report = pdbiox::validate::b_factor_distribution(structure, &selection, threshold)
        .map_err(|error| {
            eprintln!("B-factor validation failed: {error}");
            Exit::Consistency
        })?;
    Ok(report
        .outliers
        .into_iter()
        .map(|outlier| {
            Row::new(
                "b-factor",
                outlier.atom.get().to_string(),
                format!("value={} z_score={}", outlier.value, outlier.z_score),
            )
        })
        .collect())
}

fn altloc_rows(
    structure: &pdbiox::Structure,
    options: &ValidationOptions<'_>,
    context: Context,
) -> Result<Vec<Row>, Exit> {
    let (Some(expected_sum), Some(tolerance)) =
        (options.altloc_expected_sum, options.altloc_tolerance)
    else {
        eprintln!("altloc validation requires --altloc-expected-sum and --altloc-tolerance");
        return Err(Exit::Usage);
    };
    let report = pdbiox::validate::altloc_occupancy_sums(
        structure,
        context.policy.identifiers,
        pdbiox::validate::AltlocOccupancyOptions {
            expected_sum,
            tolerance,
        },
    )
    .map_err(|error| {
        eprintln!("alternate occupancy validation failed: {error}");
        Exit::Consistency
    })?;
    Ok(report
        .records
        .into_iter()
        .filter_map(|record| {
            record.issue.map(|issue| {
                Row::new(
                    "altloc-occupancy",
                    format!("{}:{}", record.residue.get(), record.atom_name),
                    format!("{issue:?}"),
                )
            })
        })
        .collect())
}

fn ccd_rows(
    structure: &pdbiox::Structure,
    options: &ValidationOptions<'_>,
    context: Context,
) -> Result<Vec<Row>, Exit> {
    let source = crate::chemistry::resolve_ccd(options.ccd, options.ccd_version, context)?;
    let provider = crate::chemistry::load_ccd(source.path, source.version, source.context)?;
    let report = pdbiox::validate::ccd_missing_atoms(structure, &provider, source.context.policy)
        .map_err(|finding| {
        source.context.findings(&[finding], "ccd-completeness");
        Exit::Consistency
    })?;
    Ok(report
        .residues
        .into_iter()
        .filter(|residue| !residue.missing.is_empty() || !residue.ambiguous.is_empty())
        .map(|residue| {
            Row::new(
                "ccd-completeness",
                residue.residue.get().to_string(),
                format!(
                    "missing={} ambiguous={}",
                    residue.missing.join(","),
                    residue.ambiguous.join(",")
                ),
            )
        })
        .collect())
}

fn core_rows(structure: &pdbiox::Structure, context: Context, input: &Path) -> Vec<Row> {
    let findings = pdbiox::validate(structure.data());
    context.findings(&findings, &input.display().to_string());
    findings
        .iter()
        .map(|finding| Row::new("core", finding.code().to_string(), "violation"))
        .collect()
}

fn quality_rows(structure: &pdbiox::Structure) -> Vec<Row> {
    pdbiox::validate::quality_flags(structure)
        .into_iter()
        .map(|flag| {
            Row::new(
                "quality",
                flag.atom.to_string(),
                format!("{:?}", flag.issue),
            )
        })
        .collect()
}

fn geometry_rows(
    structure: &pdbiox::Structure,
    options: &ValidationOptions<'_>,
) -> Result<Vec<Row>, Exit> {
    let (
        Some(bond_tolerance),
        Some(planarity_tolerance),
        Some(plane_relative_tolerance),
        Some(plane_maximum_sweeps),
    ) = (
        options.bond_tolerance,
        options.planarity_tolerance,
        options.plane_relative_tolerance,
        options.plane_maximum_sweeps,
    )
    else {
        eprintln!(
            "geometry validation requires --bond-tolerance, --planarity-tolerance, \
             --plane-relative-tolerance and --plane-maximum-sweeps"
        );
        return Err(Exit::Usage);
    };
    let mut rows: Vec<Row> = pdbiox::validate::bond_length_deviations(structure, bond_tolerance)
        .into_iter()
        .map(|flag| {
            Row::new(
                "bond-length",
                format!("{}-{}", flag.atom_a, flag.atom_b),
                flag.deviation.to_string(),
            )
        })
        .collect();
    let planarity = pdbiox::validate::nonplanar_aromatic_rings(
        structure,
        pdbiox::validate::PlanarityOptions {
            maximum_deviation: planarity_tolerance,
            plane_fit: pdbiox::EigenOptions {
                relative_tolerance: plane_relative_tolerance,
                maximum_sweeps: plane_maximum_sweeps,
            },
        },
    )
    .map_err(|error| {
        eprintln!("aromatic planarity validation failed: {error}");
        Exit::Consistency
    })?;
    rows.extend(planarity.into_iter().map(|flag| {
        Row::new(
            "planarity",
            flag.residue.to_string(),
            flag.deviation.to_string(),
        )
    }));
    Ok(rows)
}

fn clash_rows(
    structure: &pdbiox::Structure,
    options: &ValidationOptions<'_>,
) -> Result<Vec<Row>, Exit> {
    let (Some(tolerance), Some(radii)) = (options.clash_tolerance, options.radii) else {
        eprintln!("clash validation requires --clash-tolerance and --radii");
        return Err(Exit::Usage);
    };
    match pdbiox::validate::clashes(structure, tolerance, radii, SpatialBackend::Auto) {
        Ok(flags) => Ok(flags
            .into_iter()
            .map(|flag| {
                Row::new(
                    "clash",
                    format!("{}-{}", flag.first, flag.second),
                    flag.overlap.to_string(),
                )
            })
            .collect()),
        Err(error) => {
            eprintln!("clash validation failed: {error}");
            Err(Exit::Consistency)
        }
    }
}

fn completeness_rows(structure: &pdbiox::Structure, context: Context) -> Result<Vec<Row>, Exit> {
    let report = match pdbiox::validate::completeness(structure, context.policy.identifiers) {
        Ok(report) => report,
        Err(error) => {
            eprintln!("completeness validation failed: {error}");
            return Err(Exit::Consistency);
        }
    };
    Ok(report
        .into_iter()
        .filter(|chain| !chain.missing.is_empty())
        .map(|chain| {
            Row::new(
                "completeness",
                chain.chain,
                format!("{} missing of {}", chain.missing.len(), chain.canonical),
            )
        })
        .collect())
}

#[derive(Debug)]
struct Row {
    check: String,
    item: String,
    value: String,
}

impl Row {
    fn new(check: impl Into<String>, item: impl Into<String>, value: impl Into<String>) -> Self {
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
        let delimiter = match context.delimiter() {
            Some(value) => value,
            None => '\t',
        };
        let mut table = Table::new(delimiter, &["check", "item", "value"]);
        for row in rows {
            table.row([row.check.as_str(), row.item.as_str(), row.value.as_str()]);
        }
        context.result(&table.finish());
    }
}
