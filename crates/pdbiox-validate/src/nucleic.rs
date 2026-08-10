//! Per-residue nucleic geometry from CCD components and explicit atom roles.

use pdbiox_analysis::{
    NucleicTorsionError, NucleicTorsions, Pucker, nucleic_torsions, sugar_pucker,
};
use pdbiox_chem::{ComponentProvider, PolymerAtomRole, PolymerRoleProfile};
use pdbiox_core::index::ResidueIndex;
use pdbiox_core::structure::{AtomRef, ResidueRef, Structure};
use pdbiox_core::{AtomAnnotation, Diagnostic, Presence};

/// Explicit thresholds for nucleotide geometry validation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NucleicGeometryPolicy {
    /// Maximum base-atom RMS departure from a best-fit plane.
    pub maximum_base_plane_deviation: f64,
    /// Inclusive glycosidic bond-length interval.
    pub glycosidic_bond_range: [f64; 2],
    /// Inclusive phosphodiester-link interval.
    pub phosphodiester_bond_range: [f64; 2],
    /// Numerical controls for fitting nucleotide base planes.
    pub plane_fit: pdbiox_geom::EigenOptions,
}

/// One actionable nucleotide-geometry finding.
#[derive(Clone, Debug, PartialEq)]
pub enum NucleicGeometryIssue {
    /// Observed base-role atoms do not cover the explicit CCD/profile expectation.
    MissingBaseAtoms {
        /// Available base-role atoms.
        available: usize,
        /// Expected non-leaving base-role atoms.
        expected: usize,
    },
    /// Not all five explicit sugar roles were available.
    MissingSugarAtoms {
        /// Available sugar-ring roles.
        available: usize,
    },
    /// Base atoms depart too far from a plane.
    NonPlanarBase {
        /// Measured RMS departure.
        deviation: f64,
    },
    /// Glycosidic bond lies outside the configured interval.
    GlycosidicBondLength {
        /// Measured length.
        distance: f64,
    },
    /// Incoming bonded O3-to-phosphate link lies outside the configured interval.
    PhosphodiesterBondLength {
        /// Measured length.
        distance: f64,
    },
}

/// Complete measured geometry for one nucleotide residue.
#[derive(Clone, Debug, PartialEq)]
pub struct NucleicGeometryRecord {
    /// Residue being described.
    pub residue: ResidueIndex,
    /// Backbone and glycosidic torsions.
    pub torsions: NucleicTorsions,
    /// Five-membered sugar pseudorotation when all roles are present.
    pub pucker: Option<Pucker>,
    /// RMS base departure from a best-fit plane.
    pub base_plane_deviation: Option<f64>,
    /// Explicit C1-role to glycosidic-role distance.
    pub glycosidic_bond_length: Option<f64>,
    /// Bonded previous O3-role to current phosphate-role distance.
    pub incoming_phosphodiester_length: Option<f64>,
    /// Per-item validation findings.
    pub issues: Vec<NucleicGeometryIssue>,
}

/// Why nucleic geometry validation could not be evaluated.
#[derive(Debug, thiserror::Error)]
pub enum NucleicGeometryError {
    /// Thresholds are non-finite, non-positive, or reversed.
    #[error("nucleic geometry thresholds are non-finite, non-positive, or reversed")]
    InvalidPolicy,
    /// Semantic torsion projection failed.
    #[error(transparent)]
    Torsions(#[from] NucleicTorsionError),
    /// A base plane could not be decomposed.
    #[error("nucleic base plane fitting failed: {0:?}")]
    Geometry(pdbiox_geom::EigenError),
    /// Component provider access failed.
    #[error("component provider failed: {0}")]
    Provider(Diagnostic),
    /// A nucleotide residue has no deposited component identifier.
    #[error("residue {residue} has no component identity for CCD lookup")]
    MissingComponentIdentity {
        /// Affected residue.
        residue: ResidueIndex,
    },
    /// A role-annotated nucleotide lacks its explicit CCD component.
    #[error("CCD component {component} for residue {residue} is unavailable")]
    MissingComponent {
        /// Affected residue.
        residue: ResidueIndex,
        /// Deposited component ID.
        component: Box<str>,
    },
    /// A role expected to be unique occurs more than once.
    #[error("residue {residue} has multiple atoms for polymer role {role}")]
    AmbiguousRole {
        /// Affected residue.
        residue: ResidueIndex,
        /// Stable role code.
        role: i64,
    },
}

/// Measures and validates nucleotides against explicit CCD/profile semantics.
///
/// # Errors
///
/// Returns [`NucleicGeometryError`] for invalid policy, unavailable CCD data,
/// absent/ambiguous roles, or torsion projection failure.
pub fn nucleic_acid_geometry(
    structure: &Structure,
    provider: &dyn ComponentProvider,
    roles: &PolymerRoleProfile,
    policy: NucleicGeometryPolicy,
) -> Result<Vec<NucleicGeometryRecord>, NucleicGeometryError> {
    validate_policy(policy)?;
    let torsions = nucleic_torsions(structure)?;
    let mut output = Vec::with_capacity(torsions.len());
    for torsions in torsions {
        let Some(residue) = structure.data().residue(torsions.residue) else {
            continue;
        };
        let component_id =
            residue
                .name()
                .ok_or(NucleicGeometryError::MissingComponentIdentity {
                    residue: residue.index(),
                })?;
        let component = provider
            .get(component_id)
            .map_err(NucleicGeometryError::Provider)?
            .ok_or_else(|| NucleicGeometryError::MissingComponent {
                residue: residue.index(),
                component: component_id.into(),
            })?;
        output.push(measure_residue(
            structure, residue, &component, roles, torsions, policy,
        )?);
    }
    Ok(output)
}

fn measure_residue(
    structure: &Structure,
    residue: ResidueRef<'_>,
    component: &pdbiox_chem::Component,
    profile: &PolymerRoleProfile,
    torsions: NucleicTorsions,
    policy: NucleicGeometryPolicy,
) -> Result<NucleicGeometryRecord, NucleicGeometryError> {
    let base_points: Vec<_> = residue
        .atoms()
        .filter(|atom| {
            atom_role(structure, *atom)
                .is_some_and(|role| role.intersects(PolymerAtomRole::NUCLEIC_BASE_GROUP))
        })
        .filter_map(AtomRef::position)
        .collect();
    let expected_base_atoms = component
        .atoms
        .iter()
        .filter(|atom| {
            !atom.leaving
                && profile
                    .role_for(component, &atom.name)
                    .intersects(PolymerAtomRole::NUCLEIC_BASE_GROUP)
        })
        .count();
    let sugar = sugar_positions(structure, residue)?;
    let sugar_points: Vec<_> = sugar.into_iter().flatten().collect();
    let pucker = sugar_pseudorotation(sugar);
    let base_plane_deviation =
        pdbiox_geom::plane_deviation_with_options(&base_points, policy.plane_fit)
            .map_err(NucleicGeometryError::Geometry)?;
    let glycosidic_bond_length = role_position(structure, residue, PolymerAtomRole::NUCLEIC_C1)?
        .zip(role_position(
            structure,
            residue,
            PolymerAtomRole::NUCLEIC_GLYCOSIDIC,
        )?)
        .map(|(first, second)| pdbiox_geom::distance(first, second));
    let incoming_phosphodiester_length = incoming_link_length(structure, residue)?;
    let mut issues = Vec::new();
    if base_points.len() != expected_base_atoms {
        issues.push(NucleicGeometryIssue::MissingBaseAtoms {
            available: base_points.len(),
            expected: expected_base_atoms,
        });
    }
    if sugar_points.len() != sugar.len() {
        issues.push(NucleicGeometryIssue::MissingSugarAtoms {
            available: sugar_points.len(),
        });
    }
    if let Some(deviation) =
        base_plane_deviation.filter(|value| *value > policy.maximum_base_plane_deviation)
    {
        issues.push(NucleicGeometryIssue::NonPlanarBase { deviation });
    }
    if let Some(distance) =
        glycosidic_bond_length.filter(|value| !inside(*value, policy.glycosidic_bond_range))
    {
        issues.push(NucleicGeometryIssue::GlycosidicBondLength { distance });
    }
    if let Some(distance) = incoming_phosphodiester_length
        .filter(|value| !inside(*value, policy.phosphodiester_bond_range))
    {
        issues.push(NucleicGeometryIssue::PhosphodiesterBondLength { distance });
    }
    Ok(NucleicGeometryRecord {
        residue: residue.index(),
        torsions,
        pucker,
        base_plane_deviation,
        glycosidic_bond_length,
        incoming_phosphodiester_length,
        issues,
    })
}

fn sugar_positions(
    structure: &Structure,
    residue: ResidueRef<'_>,
) -> Result<[Option<[f32; 3]>; 5], NucleicGeometryError> {
    Ok([
        role_position(structure, residue, PolymerAtomRole::NUCLEIC_O4)?,
        role_position(structure, residue, PolymerAtomRole::NUCLEIC_C1)?,
        role_position(structure, residue, PolymerAtomRole::NUCLEIC_C2)?,
        role_position(structure, residue, PolymerAtomRole::NUCLEIC_C3)?,
        role_position(structure, residue, PolymerAtomRole::NUCLEIC_C4)?,
    ])
}

fn sugar_pseudorotation(points: [Option<[f32; 3]>; 5]) -> Option<Pucker> {
    let [Some(o4), Some(c1), Some(c2), Some(c3), Some(c4)] = points else {
        return None;
    };
    let torsions = [
        pdbiox_geom::dihedral(o4, c1, c2, c3),
        pdbiox_geom::dihedral(c1, c2, c3, c4),
        pdbiox_geom::dihedral(c2, c3, c4, o4),
        pdbiox_geom::dihedral(c3, c4, o4, c1),
        pdbiox_geom::dihedral(c4, o4, c1, c2),
    ];
    let values: Option<Vec<_>> = torsions
        .into_iter()
        .map(|value| value.map(pdbiox_geom::degrees))
        .collect();
    Some(sugar_pucker(values?.try_into().ok()?))
}

fn incoming_link_length(
    structure: &Structure,
    residue: ResidueRef<'_>,
) -> Result<Option<f64>, NucleicGeometryError> {
    let Some(phosphate) = role_atom(structure, residue, PolymerAtomRole::NUCLEIC_PHOSPHATE)? else {
        return Ok(None);
    };
    if !structure.data().bonds.is_available() {
        return Ok(None);
    }
    let adjacency = structure.data().bonds.adjacency(structure.atom_count());
    for neighbour in adjacency.neighbours(phosphate.index()) {
        let Some(atom) = structure.atom(*neighbour) else {
            continue;
        };
        if atom
            .residue()
            .is_some_and(|other| other.index() != residue.index())
            && atom_role(structure, atom)
                .is_some_and(|role| role.intersects(PolymerAtomRole::NUCLEIC_O3))
        {
            return Ok(atom
                .position()
                .zip(phosphate.position())
                .map(|(first, second)| pdbiox_geom::distance(first, second)));
        }
    }
    Ok(None)
}

fn role_position(
    structure: &Structure,
    residue: ResidueRef<'_>,
    required: PolymerAtomRole,
) -> Result<Option<[f32; 3]>, NucleicGeometryError> {
    Ok(role_atom(structure, residue, required)?.and_then(AtomRef::position))
}

fn role_atom<'a>(
    structure: &'a Structure,
    residue: ResidueRef<'a>,
    required: PolymerAtomRole,
) -> Result<Option<AtomRef<'a>>, NucleicGeometryError> {
    let mut matching = residue
        .atoms()
        .filter(|atom| atom_role(structure, *atom).is_some_and(|role| role.intersects(required)));
    let first = matching.next();
    if matching.next().is_some() {
        Err(NucleicGeometryError::AmbiguousRole {
            residue: residue.index(),
            role: required.code(),
        })
    } else {
        Ok(first)
    }
}

fn atom_role(structure: &Structure, atom: AtomRef<'_>) -> Option<PolymerAtomRole> {
    let AtomAnnotation::Integer(column) = structure
        .annotations()
        .get(pdbiox_core::POLYMER_ATOM_ROLE_ANNOTATION)?
    else {
        return None;
    };
    column
        .get(atom.index().get())
        .filter(|(_, presence)| *presence == Presence::Present)
        .and_then(|(code, _)| PolymerAtomRole::from_code(code))
}

fn inside(value: f64, range: [f64; 2]) -> bool {
    (range[0]..=range[1]).contains(&value)
}

fn validate_policy(policy: NucleicGeometryPolicy) -> Result<(), NucleicGeometryError> {
    let ranges = [
        policy.glycosidic_bond_range,
        policy.phosphodiester_bond_range,
    ];
    if !policy.maximum_base_plane_deviation.is_finite()
        || policy.maximum_base_plane_deviation < 0.0
        || ranges.iter().any(|range| {
            !range[0].is_finite() || !range[1].is_finite() || range[0] <= 0.0 || range[1] < range[0]
        })
    {
        Err(NucleicGeometryError::InvalidPolicy)
    } else {
        Ok(())
    }
}

#[cfg(test)]
#[path = "nucleic_tests.rs"]
mod tests;
