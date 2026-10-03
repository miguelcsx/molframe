//! Checks over bonding, stereochemistry and polymer geometry.
//!
//! The ones that name backbone or nucleotide atoms do so through the caller's
//! polymer role profile, never through atom names.

use super::{Row, Subject};
use crate::exit::Exit;
use std::path::Path;

/// The structure's bond table, or a usage error saying how to obtain one.
///
/// The kernels return no findings when bonds are absent, which would read as
/// a pass; so absence is refused here instead.
fn require_bonds(subject: &Subject<'_>, check: &str) -> Result<(), Exit> {
    if subject.structure.engine().data().bonds.is_available() {
        return Ok(());
    }
    eprintln!(
        "{check} validation needs a bond table; read a format that carries connectivity or \
         validate with a chemistry-annotated structure"
    );
    Err(Exit::Usage)
}

fn required<T>(value: Option<T>, check: &str, flags: &str) -> Result<T, Exit> {
    value.ok_or_else(|| {
        eprintln!("{check} validation requires {flags}");
        Exit::Usage
    })
}

pub(super) fn valence_rows(subject: &Subject<'_>) -> Result<Vec<Row>, Exit> {
    require_bonds(subject, "valence")?;
    Ok(
        molframe::validation::overvalent_atoms(subject.structure.engine())
            .into_iter()
            .map(|error| {
                Row::new(
                    "valence",
                    error.atom.get().to_string(),
                    format!("bonds={} maximum={}", error.bonds, error.maximum),
                )
            })
            .collect(),
    )
}

pub(super) fn ligand_rows(subject: &Subject<'_>) -> Result<Vec<Row>, Exit> {
    require_bonds(subject, "ligand")?;
    let tolerance = required(
        subject.options.polymer.ligand_bond_tolerance,
        "ligand",
        "--ligand-bond-tolerance",
    )?;
    let report = molframe::validation::ligand_geometry(subject.structure.engine(), tolerance);
    eprintln!(
        "ligand geometry: assessed {} of {} hetero bonds",
        report.assessed, report.intended
    );
    Ok(report
        .outliers
        .into_iter()
        .map(|flag| {
            Row::new(
                "ligand-bond",
                format!("{}-{}", flag.atom_a, flag.atom_b),
                flag.deviation.to_string(),
            )
        })
        .collect())
}

pub(super) fn chirality_rows(subject: &Subject<'_>) -> Result<Vec<Row>, Exit> {
    let minimum_abs_volume = required(
        subject.options.polymer.chirality_minimum_volume,
        "chirality",
        "--chirality-minimum-volume",
    )?;
    let provider = provider(subject)?;
    let report = molframe::validation::chirality_outliers(
        subject.structure.engine(),
        &provider,
        subject.context.policy,
        molframe::validation::ChiralityOptions { minimum_abs_volume },
    )
    .map_err(|finding| {
        subject.context.findings(&[finding], "chirality");
        Exit::Consistency
    })?;
    subject.context.findings(&report.findings, "chirality");
    eprintln!(
        "chirality: assessed {} of {} stereocentres (CCD {})",
        report.assessed,
        report.intended,
        report.dictionary_version.as_str()
    );
    Ok(report
        .flags
        .into_iter()
        .map(|flag| {
            Row::new(
                "chirality",
                format!("{}:{}", flag.residue.get(), flag.centre.get()),
                format!(
                    "expected={:?} issue={:?} observed_volume={}",
                    flag.expected, flag.issue, flag.observed_volume
                ),
            )
        })
        .collect())
}

pub(super) fn cis_peptide_rows(subject: &Subject<'_>) -> Result<Vec<Row>, Exit> {
    let threshold = required(
        subject.options.polymer.cis_threshold_degrees,
        "cis-peptide",
        "--cis-threshold-degrees",
    )?;
    let roles = with_roles(subject, "cis-peptide")?;
    let found = molframe::validation::cis_peptides(&roles, threshold).map_err(|finding| {
        subject.context.findings(&[finding], "cis-peptide");
        Exit::Consistency
    })?;
    Ok(found
        .into_iter()
        .map(|peptide| {
            Row::new(
                "cis-peptide",
                peptide.residue.get().to_string(),
                peptide.omega.to_string(),
            )
        })
        .collect())
}

pub(super) fn nucleic_rows(subject: &Subject<'_>) -> Result<Vec<Row>, Exit> {
    use molframe::validation::{NucleicGeometryIssue as Issue, NucleicGeometryPolicy};
    let polymer = &subject.options.polymer;
    let (Some(deviation), Some(glycosidic), Some(phosphodiester)) = (
        polymer.base_plane_deviation,
        polymer.glycosidic_range.as_deref(),
        polymer.phosphodiester_range.as_deref(),
    ) else {
        eprintln!(
            "nucleic validation requires --base-plane-deviation, --glycosidic-range, \
             --phosphodiester-range, --plane-relative-tolerance and --plane-maximum-sweeps"
        );
        return Err(Exit::Usage);
    };
    let (Some(relative_tolerance), Some(maximum_sweeps)) = (
        subject.options.plane_relative_tolerance,
        subject.options.plane_maximum_sweeps,
    ) else {
        eprintln!(
            "nucleic validation requires --plane-relative-tolerance and --plane-maximum-sweeps"
        );
        return Err(Exit::Usage);
    };
    let (&[g_low, g_high], &[p_low, p_high]) = (glycosidic, phosphodiester) else {
        eprintln!("--glycosidic-range and --phosphodiester-range take exactly MIN MAX");
        return Err(Exit::Usage);
    };
    let profile_path = required(
        subject.options.polymer.role_profile.as_deref(),
        "nucleic",
        "--role-profile",
    )?;
    let profile = crate::role_profile::load(profile_path)?;
    let provider = provider(subject)?;
    let roles = crate::role_profile::apply(
        subject.structure,
        &profile,
        profile_path,
        &provider,
        subject.context,
    )?;
    let records = molframe::validation::nucleic_acid_geometry(
        &roles,
        &provider,
        &profile,
        NucleicGeometryPolicy {
            maximum_base_plane_deviation: deviation,
            glycosidic_bond_range: [g_low, g_high],
            phosphodiester_bond_range: [p_low, p_high],
            plane_fit: molframe::geometry::EigenOptions {
                relative_tolerance,
                maximum_sweeps,
            },
        },
    )
    .map_err(|error| {
        eprintln!("nucleic geometry validation failed: {error}");
        Exit::Consistency
    })?;
    Ok(records
        .into_iter()
        .flat_map(|record| {
            let residue = record.residue.get().to_string();
            record.issues.into_iter().map(move |issue| {
                let detail = match issue {
                    Issue::MissingBaseAtoms {
                        available,
                        expected,
                    } => format!("missing-base-atoms available={available} expected={expected}"),
                    Issue::MissingSugarAtoms { available } => {
                        format!("missing-sugar-atoms available={available}")
                    }
                    Issue::NonPlanarBase { deviation } => format!("non-planar-base {deviation}"),
                    Issue::GlycosidicBondLength { distance } => {
                        format!("glycosidic-bond-length {distance}")
                    }
                    Issue::PhosphodiesterBondLength { distance } => {
                        format!("phosphodiester-bond-length {distance}")
                    }
                };
                Row::new("nucleic", residue.clone(), detail)
            })
        })
        .collect())
}

/// The explicit dictionary shared by every dictionary-reading check.
pub(super) fn provider(subject: &Subject<'_>) -> Result<molframe::chemistry::CifProvider, Exit> {
    let (Some(path), Some(version)) = (subject.context.ccd, subject.context.ccd_version) else {
        eprintln!("this check requires --ccd and --ccd-version");
        return Err(Exit::Usage);
    };
    crate::chemistry::load_ccd(path, version, subject.context)
}

/// The structure with the caller's polymer roles applied.
pub(super) fn with_roles(
    subject: &Subject<'_>,
    check: &str,
) -> Result<molframe_core::Structure, Exit> {
    let path: &Path = required(
        subject.options.polymer.role_profile.as_deref(),
        check,
        "--role-profile",
    )?;
    let profile = crate::role_profile::load(path)?;
    let provider = provider(subject)?;
    crate::role_profile::apply(
        subject.structure,
        &profile,
        path,
        &provider,
        subject.context,
    )
}
