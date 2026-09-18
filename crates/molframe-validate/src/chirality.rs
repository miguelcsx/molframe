//! CCD-declared stereocentres compared with their reference geometry.
//!
//! Atom identity, connectivity, expected configuration and ideal coordinates
//! come from the selected component provider. The same ordered neighbours are
//! used for the reference and observed scalar triple products, so an inversion
//! changes their relative sign without implementing a second stereochemical
//! priority system. Missing chemistry is reported and never replaced by local
//! residue or atom-name conventions.

use molframe_chem::{Component, ComponentProvider, StereoConfiguration};
use molframe_core::contract::DictionaryVersion;
use molframe_core::index::{AtomIndex, ResidueIndex};
use molframe_core::structure::{ResidueRef, Structure};
use molframe_core::{AnalysisPolicy, Code, Diagnostic};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

/// Explicit numerical policy for stereocentre validation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChiralityOptions {
    /// Minimum absolute scalar triple product considered non-degenerate.
    pub minimum_abs_volume: f64,
}

/// Why a CCD-declared stereocentre failed validation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ChiralityIssue {
    /// Observed and CCD-reference geometries have opposite handedness.
    Inverted,
    /// The observed centre is too close to planar for a stable configuration.
    Degenerate,
}

/// One actionable stereochemical disagreement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChiralityFlag {
    /// Residue carrying the stereocentre.
    pub residue: ResidueIndex,
    /// Model atom at the centre.
    pub centre: AtomIndex,
    /// Configuration declared by the component dictionary.
    pub expected: StereoConfiguration,
    /// Kind of geometric disagreement.
    pub issue: ChiralityIssue,
    /// Signed volume measured in model coordinates.
    pub observed_volume: f64,
    /// Signed volume from the same ordered neighbours in CCD ideal coordinates.
    pub reference_volume: f64,
}

/// Stereochemical findings with exact chemistry and selection provenance.
#[derive(Clone, Debug)]
pub struct ChiralityReport {
    /// Sorted per-centre disagreements.
    pub flags: Vec<ChiralityFlag>,
    /// Unknown components, missing reference geometry and altloc warnings.
    pub findings: Vec<Diagnostic>,
    /// Exact component dictionary used.
    pub dictionary_version: DictionaryVersion,
    /// Numerical policy used for the decision.
    pub options: ChiralityOptions,
    /// CCD stereocentres plus unresolved component validation targets.
    pub intended: usize,
    /// Targets for which a geometric handedness decision was possible.
    pub assessed: usize,
}

/// Compares every assessable CCD-declared stereocentre with CCD ideal geometry.
///
/// Runs in `O(atoms + component bonds)` time after provider lookup caching.
/// No amino-acid-specific or atom-name fallback is performed.
///
/// # Errors
///
/// Returns a diagnostic for invalid options or component-provider failure.
pub fn chirality_outliers(
    structure: &Structure,
    provider: &dyn ComponentProvider,
    policy: &AnalysisPolicy,
    options: ChiralityOptions,
) -> Result<ChiralityReport, Diagnostic> {
    validate_options(options)?;
    let selected = structure.resolve_altlocs(policy);
    let mut findings = selected.warnings;
    let mut flags = Vec::new();
    let mut components: BTreeMap<Box<str>, Arc<Component>> = BTreeMap::new();
    let mut unresolved = BTreeSet::new();
    let mut intended = 0usize;
    let mut assessed = 0usize;
    for residue in structure.data().residues() {
        let Some(component_id) = residue.name() else {
            continue;
        };
        let component = match components.get(component_id) {
            Some(component) => Some(component.clone()),
            None => provider.get(component_id)?,
        };
        let Some(component) = component else {
            unresolved.insert(Box::<str>::from(component_id));
            intended += 1;
            continue;
        };
        components
            .entry(Box::<str>::from(component_id))
            .or_insert_with(|| component.clone());
        let coverage = assess_residue(
            residue,
            &component,
            &selected.value,
            options,
            &mut flags,
            &mut findings,
        );
        intended += coverage.0;
        assessed += coverage.1;
    }
    findings.extend(
        unresolved
            .into_iter()
            .map(|component| Diagnostic::new(Code::W3201).with_context("component", component)),
    );
    flags.sort_by_key(|flag| (flag.residue.get(), flag.centre.get()));
    Ok(ChiralityReport {
        flags,
        findings,
        dictionary_version: provider.version().clone(),
        options,
        intended,
        assessed,
    })
}

fn assess_residue(
    residue: ResidueRef<'_>,
    component: &Component,
    selected: &molframe_core::AtomSelection,
    options: ChiralityOptions,
    flags: &mut Vec<ChiralityFlag>,
    findings: &mut Vec<Diagnostic>,
) -> (usize, usize) {
    let intended = component
        .atoms
        .iter()
        .filter(|atom| {
            atom.stereo
                .is_some_and(|value| value != StereoConfiguration::Mixed)
        })
        .count();
    let Some(ideal) = component.ideal_coordinates.as_deref() else {
        report_unassessed(residue, findings, "CCD ideal coordinates are absent");
        return (intended, 0);
    };
    if ideal.len() != component.atoms.len() {
        report_unassessed(residue, findings, "CCD ideal coordinates are misaligned");
        return (intended, 0);
    }
    let observed = observed_atoms(residue, selected);
    let mut assessed = 0usize;
    for (centre, atom) in component.atoms.iter().enumerate() {
        let Some(expected) = atom
            .stereo
            .filter(|value| *value != StereoConfiguration::Mixed)
        else {
            continue;
        };
        let Some(neighbours) = assessable_neighbours(component, centre, &observed) else {
            continue;
        };
        let Some((centre_atom, centre_position)) = observed.get(atom.name.as_ref()).copied() else {
            continue;
        };
        let reference_volume = volume(ideal[centre], neighbours.map(|index| ideal[index]));
        if reference_volume.abs() < options.minimum_abs_volume {
            report_unassessed(residue, findings, "CCD stereocentre geometry is degenerate");
            continue;
        }
        let model_points = neighbours.map(|index| observed[component.atoms[index].name.as_ref()].1);
        let observed_volume = volume(centre_position, model_points);
        assessed += 1;
        let issue = if observed_volume.abs() < options.minimum_abs_volume {
            Some(ChiralityIssue::Degenerate)
        } else if observed_volume.is_sign_positive() != reference_volume.is_sign_positive() {
            Some(ChiralityIssue::Inverted)
        } else {
            None
        };
        if let Some(issue) = issue {
            flags.push(ChiralityFlag {
                residue: residue.index(),
                centre: centre_atom,
                expected,
                issue,
                observed_volume,
                reference_volume,
            });
        }
    }
    (intended, assessed)
}

fn observed_atoms<'a>(
    residue: ResidueRef<'a>,
    selected: &molframe_core::AtomSelection,
) -> BTreeMap<&'a str, (AtomIndex, [f32; 3])> {
    residue
        .atoms()
        .filter(|atom| selected.contains(atom.index().get()))
        .filter_map(|atom| Some((atom.name()?, (atom.index(), atom.position()?))))
        .collect()
}

fn assessable_neighbours(
    component: &Component,
    centre: usize,
    observed: &BTreeMap<&str, (AtomIndex, [f32; 3])>,
) -> Option<[usize; 3]> {
    let centre_name = component.atoms.get(centre)?.name.as_ref();
    if !observed.contains_key(centre_name) {
        return None;
    }
    let mut neighbours = component
        .bonds
        .iter()
        .filter_map(|bond| {
            let name = if bond.atom_a.as_ref() == centre_name {
                Some(bond.atom_b.as_ref())
            } else if bond.atom_b.as_ref() == centre_name {
                Some(bond.atom_a.as_ref())
            } else {
                None
            }?;
            component
                .atoms
                .iter()
                .position(|atom| atom.name.as_ref() == name && observed.contains_key(name))
        })
        .collect::<Vec<_>>();
    neighbours.sort_unstable();
    neighbours.dedup();
    let [first, second, third, ..] = neighbours.as_slice() else {
        return None;
    };
    Some([*first, *second, *third])
}

fn volume(centre: [f32; 3], neighbours: [[f32; 3]; 3]) -> f64 {
    let vectors = neighbours.map(|point| subtract(point, centre));
    let cross = [
        vectors[1][1] * vectors[2][2] - vectors[1][2] * vectors[2][1],
        vectors[1][2] * vectors[2][0] - vectors[1][0] * vectors[2][2],
        vectors[1][0] * vectors[2][1] - vectors[1][1] * vectors[2][0],
    ];
    vectors[0][0] * cross[0] + vectors[0][1] * cross[1] + vectors[0][2] * cross[2]
}

fn subtract(point: [f32; 3], origin: [f32; 3]) -> [f64; 3] {
    [
        f64::from(point[0]) - f64::from(origin[0]),
        f64::from(point[1]) - f64::from(origin[1]),
        f64::from(point[2]) - f64::from(origin[2]),
    ]
}

fn report_unassessed(residue: ResidueRef<'_>, findings: &mut Vec<Diagnostic>, reason: &str) {
    findings.push(
        Diagnostic::new(Code::W3205)
            .with_context("residue", residue.index().get().to_string())
            .with_context("reason", reason),
    );
}

fn validate_options(options: ChiralityOptions) -> Result<(), Diagnostic> {
    if options.minimum_abs_volume.is_finite() && options.minimum_abs_volume > 0.0 {
        Ok(())
    } else {
        Err(Diagnostic::new(Code::E4003)
            .with_context("required", "positive finite chirality volume threshold"))
    }
}

#[cfg(test)]
#[path = "chirality_tests.rs"]
mod tests;
