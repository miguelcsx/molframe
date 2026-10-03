//! The checks that need only explicit thresholds and, for some, a dictionary.

use super::{Row, ValidationOptions};
use crate::exit::Exit;
use crate::report::Context;
use molframe::spatial::SpatialBackend;
use std::path::Path;

pub(super) fn b_factor_rows(
    structure: &molframe::Structure,
    options: ValidationOptions<'_>,
) -> Result<Vec<Row>, Exit> {
    let Some(threshold) = options.b_factor_z_score else {
        eprintln!("B-factor validation requires --b-factor-z-score");
        return Err(Exit::Usage);
    };
    let selection = molframe_core::selection::AtomSelection::All(structure.atom_count());
    let report =
        molframe::validation::b_factor_distribution(structure.engine(), &selection, threshold)
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

pub(super) fn altloc_rows(
    structure: &molframe::Structure,
    options: ValidationOptions<'_>,
    context: Context,
) -> Result<Vec<Row>, Exit> {
    let (Some(expected_sum), Some(tolerance)) =
        (options.altloc_expected_sum, options.altloc_tolerance)
    else {
        eprintln!("altloc validation requires --altloc-expected-sum and --altloc-tolerance");
        return Err(Exit::Usage);
    };
    let report = molframe::validation::altloc_occupancy_sums(
        structure.engine(),
        context.policy.identifiers,
        molframe::validation::AltlocOccupancyOptions {
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

pub(super) fn ccd_rows(
    structure: &molframe::Structure,
    options: ValidationOptions<'_>,
    context: Context,
) -> Result<Vec<Row>, Exit> {
    let source = crate::chemistry::resolve_ccd(
        options.chemistry.ccd.as_deref(),
        options.chemistry.ccd_version.as_deref(),
        context,
    )?;
    let provider = crate::chemistry::load_ccd(source.path, source.version, source.context)?;
    let report = molframe::validation::ccd_missing_atoms(
        structure.engine(),
        &provider,
        source.context.policy,
    )
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

pub(super) fn core_rows(
    structure: &molframe::Structure,
    context: Context,
    input: &Path,
) -> Vec<Row> {
    let findings = molframe_core::structure::validate(structure.engine().data());
    context.findings(&findings, &input.display().to_string());
    findings
        .iter()
        .map(|finding| Row::new("core", finding.code().to_string(), "violation"))
        .collect()
}

pub(super) fn quality_rows(structure: &molframe::Structure) -> Vec<Row> {
    molframe::validation::quality_flags(structure.engine())
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

pub(super) fn geometry_rows(
    structure: &molframe::Structure,
    options: ValidationOptions<'_>,
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
    let mut rows: Vec<Row> =
        molframe::validation::bond_length_deviations(structure.engine(), bond_tolerance)
            .into_iter()
            .map(|flag| {
                Row::new(
                    "bond-length",
                    format!("{}-{}", flag.atom_a, flag.atom_b),
                    flag.deviation.to_string(),
                )
            })
            .collect();
    let planarity = molframe::validation::nonplanar_aromatic_rings(
        structure.engine(),
        molframe::validation::PlanarityOptions {
            maximum_deviation: planarity_tolerance,
            plane_fit: molframe::geometry::EigenOptions {
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

pub(super) fn clash_rows(
    structure: &molframe::Structure,
    options: ValidationOptions<'_>,
    execution: &molframe_core::ExecutionContext,
) -> Result<Vec<Row>, Exit> {
    let (Some(tolerance), Some(radii)) = (
        options.clash_tolerance,
        options
            .radii
            .map(Into::<molframe::chemistry::RadiusSet>::into),
    ) else {
        eprintln!("clash validation requires --clash-tolerance and --radii");
        return Err(Exit::Usage);
    };
    match molframe::validation::clashes(
        structure.engine(),
        tolerance,
        radii,
        SpatialBackend::Auto,
        execution,
    ) {
        Ok(flags) => Ok(flags
            .iter()
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

pub(super) fn completeness_rows(
    structure: &molframe::Structure,
    context: Context,
) -> Result<Vec<Row>, Exit> {
    let report =
        match molframe::validation::completeness(structure.engine(), context.policy.identifiers) {
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
