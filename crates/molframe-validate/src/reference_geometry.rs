//! Bond and angle validation against explicit CCD ideal coordinates.

use molframe_chem::{Component, ComponentProvider};
use molframe_core::contract::Namespace;
use molframe_core::index::{AtomIndex, ResidueIndex};
use molframe_core::structure::{ResidueRef, Structure};
use molframe_core::{Code, Diagnostic};
use std::collections::BTreeMap;

/// Explicit decision thresholds for CCD reference geometry.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReferenceGeometryOptions {
    /// Maximum absolute bond-length departure in ångström.
    pub maximum_bond_deviation: f64,
    /// Maximum absolute bond-angle departure in degrees.
    pub maximum_angle_deviation_degrees: f64,
}

/// One bond outside its CCD ideal-coordinate threshold.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReferenceBondFlag {
    /// Residue containing the bond.
    pub residue: ResidueIndex,
    /// First observed atom.
    pub first: AtomIndex,
    /// Second observed atom.
    pub second: AtomIndex,
    /// Observed bond length.
    pub observed: f64,
    /// CCD ideal-coordinate bond length.
    pub expected: f64,
    /// Signed `observed - expected` departure.
    pub deviation: f64,
}

/// One bonded angle outside its CCD ideal-coordinate threshold.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReferenceAngleFlag {
    /// Residue containing the angle.
    pub residue: ResidueIndex,
    /// First outer atom.
    pub first: AtomIndex,
    /// Central atom.
    pub centre: AtomIndex,
    /// Second outer atom.
    pub third: AtomIndex,
    /// Observed angle in degrees.
    pub observed_degrees: f64,
    /// CCD ideal-coordinate angle in degrees.
    pub expected_degrees: f64,
    /// Signed observed-minus-expected departure in degrees.
    pub deviation_degrees: f64,
}

/// CCD geometry findings and exact assessment coverage.
#[derive(Clone, Debug)]
pub struct ReferenceGeometryReport {
    /// Bond-length outliers.
    pub bonds: Vec<ReferenceBondFlag>,
    /// Bond-angle outliers.
    pub angles: Vec<ReferenceAngleFlag>,
    /// Missing components or ideal-coordinate records.
    pub findings: Vec<Diagnostic>,
    /// CCD bond and angle targets selected for validation.
    pub intended: usize,
    /// Targets with all required observed and ideal coordinates.
    pub assessed: usize,
    /// Thresholds used.
    pub options: ReferenceGeometryOptions,
}

/// Validates intracomponent bonds and angles against CCD ideal coordinates.
///
/// # Errors
///
/// Returns invalid options, explicit-namespace ambiguity or provider diagnostics.
pub fn reference_geometry(
    structure: &Structure,
    provider: &dyn ComponentProvider,
    namespace: Namespace,
    options: ReferenceGeometryOptions,
) -> Result<ReferenceGeometryReport, Diagnostic> {
    if !options.maximum_bond_deviation.is_finite()
        || options.maximum_bond_deviation < 0.0
        || !options.maximum_angle_deviation_degrees.is_finite()
        || options.maximum_angle_deviation_degrees < 0.0
    {
        return Err(Diagnostic::new(Code::E4003).with_context("options", "reference geometry"));
    }
    let mut report = ReferenceGeometryReport {
        bonds: Vec::new(),
        angles: Vec::new(),
        findings: Vec::new(),
        intended: 0,
        assessed: 0,
        options,
    };
    for residue in structure.data().residues() {
        assess_residue(residue, provider, namespace, &mut report)?;
    }
    report
        .bonds
        .sort_by_key(|flag| (flag.residue.get(), flag.first.get(), flag.second.get()));
    report.angles.sort_by_key(|flag| {
        (
            flag.residue.get(),
            flag.centre.get(),
            flag.first.get(),
            flag.third.get(),
        )
    });
    Ok(report)
}

fn assess_residue(
    residue: ResidueRef<'_>,
    provider: &dyn ComponentProvider,
    namespace: Namespace,
    report: &mut ReferenceGeometryReport,
) -> Result<(), Diagnostic> {
    let component_id = match namespace {
        Namespace::Label => residue.name(),
        Namespace::Auth => residue.auth_name(),
        _ => {
            return Err(Diagnostic::new(Code::E4003)
                .with_context("namespace", "explicit residue namespace required"));
        }
    };
    let Some(component_id) = component_id else {
        report.intended += 1;
        report
            .findings
            .push(Diagnostic::new(Code::W3201).with_context("required", "component identifier"));
        return Ok(());
    };
    let Some(component) = provider.get(component_id)? else {
        report.intended += 1;
        report
            .findings
            .push(Diagnostic::new(Code::W3201).with_context("component", component_id));
        return Ok(());
    };
    assess_component(residue, &component, namespace, report);
    Ok(())
}

fn assess_component(
    residue: ResidueRef<'_>,
    component: &Component,
    namespace: Namespace,
    report: &mut ReferenceGeometryReport,
) {
    let angles = component_angles(component);
    report.intended += component.bonds.len() + angles.len();
    let Some(ideal) = component
        .ideal_coordinates
        .as_deref()
        .filter(|coordinates| coordinates.len() == component.atoms.len())
    else {
        report.findings.push(
            Diagnostic::new(Code::W3201)
                .with_context("component", component.id.clone())
                .with_context("required", "aligned CCD ideal coordinates"),
        );
        return;
    };
    let observed = observed_atoms(residue, namespace);
    for bond in component.bonds.iter() {
        let (Some(&first), Some(&second)) = (
            observed.get(bond.atom_a.as_ref()),
            observed.get(bond.atom_b.as_ref()),
        ) else {
            continue;
        };
        let (Some(first_index), Some(second_index)) = (
            component_index(component, &bond.atom_a),
            component_index(component, &bond.atom_b),
        ) else {
            continue;
        };
        report.assessed += 1;
        assess_bond(
            residue.index(),
            first,
            second,
            ideal[first_index],
            ideal[second_index],
            report,
        );
    }
    for [first_index, centre_index, third_index] in angles {
        let names = [
            component.atoms[first_index].name.as_ref(),
            component.atoms[centre_index].name.as_ref(),
            component.atoms[third_index].name.as_ref(),
        ];
        let [Some(&first), Some(&centre), Some(&third)] = names.map(|name| observed.get(name))
        else {
            continue;
        };
        let (Some(observed_angle), Some(expected_angle)) = (
            molframe_geom::angle(first.1, centre.1, third.1),
            molframe_geom::angle(ideal[first_index], ideal[centre_index], ideal[third_index]),
        ) else {
            continue;
        };
        report.assessed += 1;
        let observed_degrees = molframe_geom::degrees(observed_angle);
        let expected_degrees = molframe_geom::degrees(expected_angle);
        let deviation_degrees = observed_degrees - expected_degrees;
        if deviation_degrees.abs() > report.options.maximum_angle_deviation_degrees {
            report.angles.push(ReferenceAngleFlag {
                residue: residue.index(),
                first: first.0,
                centre: centre.0,
                third: third.0,
                observed_degrees,
                expected_degrees,
                deviation_degrees,
            });
        }
    }
}

fn observed_atoms(
    residue: ResidueRef<'_>,
    namespace: Namespace,
) -> BTreeMap<&str, (AtomIndex, [f32; 3])> {
    residue
        .atoms()
        .filter_map(|atom| {
            let name = match namespace {
                Namespace::Label => atom.name(),
                Namespace::Auth => atom.auth_name(),
                _ => None,
            }?;
            Some((name, (atom.index(), atom.position()?)))
        })
        .collect()
}

fn component_index(component: &Component, name: &str) -> Option<usize> {
    component
        .atoms
        .iter()
        .position(|atom| atom.name.as_ref() == name)
}

fn component_angles(component: &Component) -> Vec<[usize; 3]> {
    let mut neighbours = vec![Vec::new(); component.atoms.len()];
    for bond in component.bonds.iter() {
        let (Some(first), Some(second)) = (
            component_index(component, &bond.atom_a),
            component_index(component, &bond.atom_b),
        ) else {
            continue;
        };
        neighbours[first].push(second);
        neighbours[second].push(first);
    }
    let mut angles = Vec::new();
    for (centre, adjacent) in neighbours.iter_mut().enumerate() {
        adjacent.sort_unstable();
        adjacent.dedup();
        for first in 0..adjacent.len() {
            for third in first + 1..adjacent.len() {
                angles.push([adjacent[first], centre, adjacent[third]]);
            }
        }
    }
    angles
}

fn assess_bond(
    residue: ResidueIndex,
    first: (AtomIndex, [f32; 3]),
    second: (AtomIndex, [f32; 3]),
    ideal_first: [f32; 3],
    ideal_second: [f32; 3],
    report: &mut ReferenceGeometryReport,
) {
    let observed = molframe_geom::distance(first.1, second.1);
    let expected = molframe_geom::distance(ideal_first, ideal_second);
    let deviation = observed - expected;
    if deviation.abs() > report.options.maximum_bond_deviation {
        report.bonds.push(ReferenceBondFlag {
            residue,
            first: first.0,
            second: second.0,
            observed,
            expected,
            deviation,
        });
    }
}

#[cfg(test)]
#[path = "reference_geometry_tests.rs"]
mod tests;
