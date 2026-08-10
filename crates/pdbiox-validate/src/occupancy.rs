//! Alternate-location occupancy-sum validation without residue-name rules.

use pdbiox_core::contract::Namespace;
use pdbiox_core::index::ResidueIndex;
use pdbiox_core::structure::{AtomRef, Structure};
use std::collections::BTreeMap;

/// Explicit expected sum and numerical tolerance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AltlocOccupancyOptions {
    /// Expected sum across alternate sites for one residue atom.
    pub expected_sum: f64,
    /// Maximum accepted absolute deviation from `expected_sum`.
    pub tolerance: f64,
}

/// Why one alternate-site atom group is reportable.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AltlocOccupancyIssue {
    /// At least one alternate site has no recorded occupancy.
    MissingOccupancy,
    /// The recorded sum differs from the explicit expected value.
    SumMismatch,
}

/// Validation result for one residue/atom-name alternate-site group.
#[derive(Clone, Debug, PartialEq)]
pub struct AltlocOccupancyRecord {
    /// Topology residue index.
    pub residue: ResidueIndex,
    /// Atom name in the requested namespace.
    pub atom_name: String,
    /// Number of labelled alternate sites.
    pub alternatives: usize,
    /// Number carrying a recorded occupancy.
    pub assessed: usize,
    /// Sum when every alternate site has a recorded occupancy.
    pub occupancy_sum: Option<f64>,
    /// The issue found, absent for a valid complete group.
    pub issue: Option<AltlocOccupancyIssue>,
}

/// Whole-structure coverage and per-group occupancy results.
#[derive(Clone, Debug, PartialEq)]
pub struct AltlocOccupancyReport {
    /// Alternate-site occupancy values intended for assessment.
    pub intended: usize,
    /// Alternate-site occupancy values actually recorded.
    pub assessed: usize,
    /// Stable residue/name ordered groups.
    pub records: Vec<AltlocOccupancyRecord>,
}

/// Invalid options or identifier semantics.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum AltlocOccupancyError {
    /// Expected sum and tolerance must be finite, with non-negative tolerance.
    #[error("alternate occupancy controls must be finite with non-negative tolerance")]
    InvalidOptions,
    /// Explicit namespace needs a concrete caller choice.
    #[error("alternate occupancy validation requires label or auth atom names")]
    ExplicitNamespace,
    /// A labelled alternate site lacks the requested atom name.
    #[error("an alternate site lacks its requested atom name")]
    MissingAtomName,
    /// This crate does not understand a newer namespace.
    #[error("unsupported atom-name namespace")]
    UnsupportedNamespace,
}

/// Validates occupancy sums for every labelled residue/atom-name group.
///
/// No component or water names are recognized here; grouping uses topology and
/// the caller-selected identifier namespace only.
///
/// # Errors
///
/// Returns invalid controls or identifier-policy errors.
pub fn altloc_occupancy_sums(
    structure: &Structure,
    namespace: Namespace,
    options: AltlocOccupancyOptions,
) -> Result<AltlocOccupancyReport, AltlocOccupancyError> {
    validate_options(options)?;
    let mut groups: BTreeMap<(ResidueIndex, String), Vec<Option<f32>>> = BTreeMap::new();
    for atom in structure
        .data()
        .atoms()
        .filter(|atom| atom.alt_label().is_some())
    {
        let residue = atom
            .residue()
            .ok_or(AltlocOccupancyError::MissingAtomName)?;
        let name = atom_name(atom, namespace)?
            .ok_or(AltlocOccupancyError::MissingAtomName)?
            .to_string();
        groups
            .entry((residue.index(), name))
            .or_default()
            .push(atom.occupancy());
    }
    let mut intended = 0usize;
    let mut assessed = 0usize;
    let records = groups
        .into_iter()
        .map(|((residue, atom_name), values)| {
            intended += values.len();
            let present: Vec<f64> = values
                .iter()
                .filter_map(|value| value.map(f64::from))
                .collect();
            assessed += present.len();
            occupancy_record(residue, atom_name, values.len(), &present, options)
        })
        .collect();
    Ok(AltlocOccupancyReport {
        intended,
        assessed,
        records,
    })
}

fn occupancy_record(
    residue: ResidueIndex,
    atom_name: String,
    alternatives: usize,
    present: &[f64],
    options: AltlocOccupancyOptions,
) -> AltlocOccupancyRecord {
    let complete = alternatives == present.len();
    let occupancy_sum = complete.then(|| present.iter().sum());
    let issue = if !complete {
        Some(AltlocOccupancyIssue::MissingOccupancy)
    } else if occupancy_sum
        .is_some_and(|sum: f64| (sum - options.expected_sum).abs() > options.tolerance)
    {
        Some(AltlocOccupancyIssue::SumMismatch)
    } else {
        None
    };
    AltlocOccupancyRecord {
        residue,
        atom_name,
        alternatives,
        assessed: present.len(),
        occupancy_sum,
        issue,
    }
}

fn validate_options(options: AltlocOccupancyOptions) -> Result<(), AltlocOccupancyError> {
    if options.expected_sum.is_finite() && options.tolerance.is_finite() && options.tolerance >= 0.0
    {
        Ok(())
    } else {
        Err(AltlocOccupancyError::InvalidOptions)
    }
}

fn atom_name(
    atom: AtomRef<'_>,
    namespace: Namespace,
) -> Result<Option<&str>, AltlocOccupancyError> {
    match namespace {
        Namespace::Label => Ok(atom.name()),
        Namespace::Auth => Ok(atom.auth_name()),
        Namespace::Explicit => Err(AltlocOccupancyError::ExplicitNamespace),
        _ => Err(AltlocOccupancyError::UnsupportedNamespace),
    }
}

#[cfg(test)]
#[path = "occupancy_tests.rs"]
mod tests;
