//! Checks that compare a structure with caller-supplied references.
//!
//! The reference libraries, rotamer profiles, plane restraints and TLS groups
//! are read by the library; this file only routes them to the kernels.

use super::structural::{provider, with_roles};
use super::{Row, Subject};
use crate::exit::Exit;
use molframe::validation::{
    PlanarityOptions, RamachandranBasin, RamachandranOptions, RotamerOptions,
};
use std::path::Path;

fn required<T>(value: Option<T>, check: &str, flags: &str) -> Result<T, Exit> {
    value.ok_or_else(|| {
        eprintln!("{check} validation requires {flags}");
        Exit::Usage
    })
}

fn input_error(what: &str, error: &dyn std::fmt::Display) -> Exit {
    eprintln!("{what} could not be read: {error}");
    Exit::Policy
}

fn library(
    subject: &Subject<'_>,
    check: &str,
) -> Result<molframe::validation::ReferenceLibrary, Exit> {
    let path: &Path = required(
        subject.options.references.reference_library.as_deref(),
        check,
        "--reference-library",
    )?;
    molframe::read_reference_library(path).map_err(|error| input_error("reference library", &error))
}

pub(super) fn ramachandran_rows(subject: &Subject<'_>) -> Result<Vec<Row>, Exit> {
    let references = &subject.options.references;
    let minimum = required(
        references.ramachandran_minimum_probability,
        "ramachandran",
        "--ramachandran-minimum-probability",
    )?;
    if references.basin.is_empty() {
        eprintln!("ramachandran validation requires at least one --basin REGION=DISTRIBUTION");
        return Err(Exit::Usage);
    }
    let library = library(subject, "ramachandran")?;
    let basins = references
        .basin
        .iter()
        .map(|basin| RamachandranBasin::new(basin.region, basin.distribution.as_str()))
        .collect::<Result<Vec<_>, _>>()
        .map_err(|error| input_error("Ramachandran basin", &error))?;
    let options = RamachandranOptions::new(&library, basins, minimum)
        .map_err(|error| input_error("Ramachandran policy", &error))?;
    let roles = with_roles(subject, "ramachandran")?;
    let outliers =
        molframe::validation::ramachandran_outliers(&roles, &options).map_err(|error| {
            eprintln!("Ramachandran validation failed: {error}");
            Exit::Consistency
        })?;
    eprintln!(
        "ramachandran: reference {} {}",
        library.id(),
        library.version()
    );
    Ok(outliers
        .into_iter()
        .map(|record| {
            Row::new(
                "ramachandran",
                record.residue.get().to_string(),
                format!(
                    "phi={} psi={} probability={}",
                    record.phi, record.psi, record.reference.probability
                ),
            )
        })
        .collect())
}

pub(super) fn rotamer_rows(subject: &Subject<'_>) -> Result<Vec<Row>, Exit> {
    let references = &subject.options.references;
    let minimum_probability = required(
        references.rotamer_minimum_probability,
        "rotamer",
        "--rotamer-minimum-probability",
    )?;
    let profile_path: &Path = required(
        references.rotamer_profile.as_deref(),
        "rotamer",
        "--rotamer-profile",
    )?;
    let profile = molframe::read_rotamer_profile(profile_path)
        .map_err(|error| input_error("rotamer profile", &error))?;
    let library = library(subject, "rotamer")?;
    let provider = provider(subject)?;
    let report = molframe::validation::rotamer_outliers(
        subject.structure.engine(),
        &provider,
        subject.context.policy,
        &library,
        &profile,
        RotamerOptions {
            minimum_probability,
        },
    )
    .map_err(|error| {
        eprintln!("rotamer validation failed: {error}");
        Exit::Consistency
    })?;
    subject.context.findings(&report.findings, "rotamer");
    eprintln!(
        "rotamer: assessed {} of {} torsions (profile {} {}, reference {} {})",
        report.assessed,
        report.intended,
        report.profile_id,
        report.profile_version,
        report.reference_id,
        report.reference_version
    );
    Ok(report
        .flags
        .into_iter()
        .map(|flag| {
            Row::new(
                "rotamer",
                format!("{}:chi{}", flag.residue.get(), flag.chi_index),
                format!(
                    "chi={} probability={} distribution={}",
                    flag.chi_degrees, flag.probability, flag.distribution
                ),
            )
        })
        .collect())
}

pub(super) fn plane_restraint_rows(subject: &Subject<'_>) -> Result<Vec<Row>, Exit> {
    use molframe::ValidationExt as _;
    let path: &Path = required(
        subject.options.references.plane_restraints.as_deref(),
        "plane-restraints",
        "--plane-restraints",
    )?;
    let (Some(maximum_deviation), Some(relative_tolerance), Some(maximum_sweeps)) = (
        subject.options.planarity_tolerance,
        subject.options.plane_relative_tolerance,
        subject.options.plane_maximum_sweeps,
    ) else {
        eprintln!(
            "plane-restraints validation requires --planarity-tolerance, \
             --plane-relative-tolerance and --plane-maximum-sweeps"
        );
        return Err(Exit::Usage);
    };
    let restraints =
        molframe::read_plane_restraints(path, subject.structure, subject.context.policy)
            .map_err(|error| input_error("plane restraints", &error))?;
    let report = subject
        .structure
        .plane_restraint_outliers(
            &restraints,
            PlanarityOptions {
                maximum_deviation,
                plane_fit: molframe::geometry::EigenOptions {
                    relative_tolerance,
                    maximum_sweeps,
                },
            },
        )
        .map_err(|error| {
            eprintln!("plane restraint validation failed: {error}");
            Exit::Consistency
        })?;
    eprintln!(
        "plane-restraints: assessed {} of {} restraints",
        report.assessed, report.intended
    );
    Ok(report
        .flags
        .into_iter()
        .map(|flag| Row::new("plane-restraint", flag.id, flag.deviation.to_string()))
        .collect())
}

pub(super) fn tls_rows(subject: &Subject<'_>) -> Result<Vec<Row>, Exit> {
    let path: &Path = required(
        subject.options.references.tls_groups.as_deref(),
        "tls",
        "--tls-groups",
    )?;
    let (Some(maximum), Some(symmetry)) = (
        subject.options.references.tls_max_deviation,
        subject.options.references.tls_symmetry_tolerance,
    ) else {
        eprintln!("tls validation requires --tls-max-deviation and --tls-symmetry-tolerance");
        return Err(Exit::Usage);
    };
    let groups = molframe::read_tls_groups(path, subject.structure, subject.context.policy)
        .map_err(|error| input_error("TLS groups", &error))?;
    let report = molframe::validation::tls_b_factor_consistency(
        subject.structure.engine(),
        &groups,
        maximum,
        symmetry,
    )
    .map_err(|error| {
        eprintln!("TLS validation failed: {error}");
        Exit::Consistency
    })?;
    eprintln!(
        "tls: assessed {} of {} atom memberships",
        report.assessed, report.intended
    );
    Ok(report
        .flags
        .into_iter()
        .map(|flag| {
            Row::new(
                "tls",
                format!("{}:{}", flag.group, flag.atom.get()),
                format!(
                    "observed={} predicted={} deviation={}",
                    flag.observed, flag.predicted, flag.deviation
                ),
            )
        })
        .collect())
}
