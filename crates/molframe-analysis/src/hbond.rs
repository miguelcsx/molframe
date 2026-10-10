//! CCD-perceived hydrogen bonds with explicit hydrogen geometry.

use molframe_core::column::Presence;
use molframe_core::index::AtomIndex;
use molframe_core::selection::AtomSelection;
use molframe_core::structure::{AtomRef, Structure};
use molframe_core::{ExecutionContext, annotation::AtomAnnotation};
use molframe_spatial::{
    PairQuery, PeriodicBox, SpatialBackend, SpatialError, SpatialSearchOptions,
    reduce_pairs_within_unsorted,
};

use crate::numeric::f64_to_f32;

/// Explicit geometric and spatial policy for hydrogen-bond detection.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HydrogenBondOptions {
    /// Maximum donor–acceptor heavy-atom distance in ångström.
    pub maximum_donor_acceptor_distance: f32,
    /// Minimum donor–hydrogen–acceptor angle in degrees.
    pub minimum_angle_degrees: f64,
    /// Spatial implementation.
    pub backend: SpatialBackend,
    /// Apply the structure unit cell and minimum-image convention.
    pub periodic: bool,
}

/// Optional refinements to hydrogen-bond detection that sit beside the core
/// geometric options.
///
/// The default preserves the behaviour of a plain [`HydrogenBondOptions`]
/// request: no hydrogen-to-acceptor cutoff and no placeholder cell accepted.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HydrogenBondPolicy {
    /// Maximum hydrogen–acceptor distance in ångström.
    ///
    /// `None` (the default) applies only the donor–acceptor distance and the
    /// angle, which keeps long, bent contacts that a hydrogen-distance cutoff
    /// would call weak. Set it (2.5 Å is a common choice) to drop them.
    pub maximum_hydrogen_acceptor_distance: Option<f64>,
    /// Accept a placeholder unit cell for periodic detection.
    ///
    /// A unit cube with right angles is what files write when they have no
    /// cell, so periodic images built from it are fictitious. Periodic
    /// requests are rejected unless this is set explicitly.
    pub allow_placeholder_cell: bool,
}

/// One fully oriented hydrogen bond.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HydrogenBond {
    /// CCD-perceived donor heavy atom.
    pub donor: AtomIndex,
    /// Explicit hydrogen covalently attached to the donor.
    pub hydrogen: AtomIndex,
    /// CCD-perceived acceptor atom.
    pub acceptor: AtomIndex,
    /// Donor–acceptor distance in ångström.
    pub donor_acceptor_distance: f32,
    /// Hydrogen–acceptor distance in ångström.
    pub hydrogen_acceptor_distance: f32,
    /// Donor–hydrogen–acceptor angle in degrees.
    pub angle_degrees: f64,
}

define_soa_table! {
    /// Native columnar storage for hydrogen-bond records.
    pub struct HydrogenBondTable for HydrogenBond {
        /// Donor atom indices.
        donor: AtomIndex,
        /// Explicit hydrogen atom indices.
        hydrogen: AtomIndex,
        /// Acceptor atom indices.
        acceptor: AtomIndex,
        /// Donor-to-acceptor distances.
        donor_acceptor_distance: f32,
        /// Hydrogen-to-acceptor distances.
        hydrogen_acceptor_distance: f32,
        /// Donor-hydrogen-acceptor angles.
        angle_degrees: f64,
    }
}

/// Hydrogen-bond detection failure.
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum HydrogenBondError {
    /// Options contain an invalid distance or angle.
    #[error("invalid hydrogen-bond geometry options")]
    InvalidOptions,
    /// CCD donor/acceptor annotations or resolved bonds are absent.
    #[error("hydrogen bonds require CCD donor/acceptor annotations and connectivity")]
    MissingChemistry,
    /// Periodic geometry was requested without a usable unit cell: the cell is
    /// absent, or it is the placeholder some files write when they have none
    /// and the policy does not allow it.
    #[error("periodic hydrogen bonds require a real unit cell (absent or placeholder cell)")]
    MissingCell,
    /// Spatial indexing rejected the request.
    #[error(transparent)]
    Spatial(#[from] SpatialError),
}

/// Finds hydrogen bonds from CCD roles, explicit hydrogens and angular geometry.
///
/// No heavy-atom-only direction or residue-name fallback is performed. Atoms
/// from incompatible alternate locations are never paired: a blank label is
/// compatible with anything and two equal labels are compatible. Results are
/// sorted by `(donor, hydrogen, acceptor)`.
///
/// # Errors
///
/// Returns [`HydrogenBondError`] for invalid options, absent chemistry,
/// unavailable periodic cell data or a spatial-search failure.
pub fn hydrogen_bonds(
    structure: &Structure,
    options: HydrogenBondOptions,
    context: &ExecutionContext,
) -> Result<HydrogenBondTable, HydrogenBondError> {
    hydrogen_bonds_with_policy(structure, options, HydrogenBondPolicy::default(), context)
}

/// [`hydrogen_bonds`] with the optional [`HydrogenBondPolicy`] refinements.
///
/// # Errors
///
/// As [`hydrogen_bonds`]; a periodic request over a placeholder cell the
/// policy does not allow is reported as [`HydrogenBondError::MissingCell`].
pub fn hydrogen_bonds_with_policy(
    structure: &Structure,
    options: HydrogenBondOptions,
    policy: HydrogenBondPolicy,
    context: &ExecutionContext,
) -> Result<HydrogenBondTable, HydrogenBondError> {
    validate_options(options, policy)?;
    if !structure.data().bonds.is_available() {
        return Err(HydrogenBondError::MissingChemistry);
    }
    let donors = role_selection(structure, molframe_core::HBOND_DONOR_ANNOTATION)?;
    let acceptors = role_selection(structure, molframe_core::HBOND_ACCEPTOR_ANNOTATION)?;
    let periodic_box = periodic_box(structure, options, policy)?;
    let adjacency = structure.data().bonds.adjacency(structure.atom_count());
    let mut output = Vec::new();

    // Only geometrically valid bonds survive, so candidate pairs are consumed as
    // they arrive instead of being collected first.
    let query = PairQuery {
        positions: structure.positions(),
        left: &donors,
        right: &acceptors,
        cutoff: options.maximum_donor_acceptor_distance,
        options: SpatialSearchOptions::with_backend(options.backend),
        periodic: periodic_box.as_ref(),
        context,
    };
    let parts =
        reduce_pairs_within_unsorted(&query, Vec::new, |output: &mut Vec<HydrogenBond>, pair| {
            for (donor, acceptor) in orientations(pair.first, pair.second, &donors, &acceptors) {
                if donor == acceptor
                    || !altlocs_compatible(structure, donor, acceptor)
                    || adjacency
                        .neighbours(AtomIndex::new(donor))
                        .binary_search(&AtomIndex::new(acceptor))
                        .is_ok()
                {
                    continue;
                }
                for hydrogen in donor_hydrogens(structure, adjacency, donor) {
                    if let Some(bond) =
                        measure(structure, donor, hydrogen, acceptor, periodic_box.as_ref())
                        && bond.angle_degrees >= options.minimum_angle_degrees
                        && policy
                            .maximum_hydrogen_acceptor_distance
                            .is_none_or(|limit| f64::from(bond.hydrogen_acceptor_distance) <= limit)
                    {
                        output.push(bond);
                    }
                }
            }
        })?;

    for part in parts {
        output.extend(part);
    }

    output.sort_by_key(|bond| (bond.donor.get(), bond.hydrogen.get(), bond.acceptor.get()));
    output.dedup_by_key(|bond| (bond.donor.get(), bond.hydrogen.get(), bond.acceptor.get()));
    Ok(output.into_iter().collect())
}

fn validate_options(
    options: HydrogenBondOptions,
    policy: HydrogenBondPolicy,
) -> Result<(), HydrogenBondError> {
    if policy
        .maximum_hydrogen_acceptor_distance
        .is_some_and(|limit| !limit.is_finite() || limit <= 0.0)
    {
        return Err(HydrogenBondError::InvalidOptions);
    }
    if options.maximum_donor_acceptor_distance.is_finite()
        && options.maximum_donor_acceptor_distance > 0.0
        && options.minimum_angle_degrees.is_finite()
        && (0.0..=180.0).contains(&options.minimum_angle_degrees)
    {
        Ok(())
    } else {
        Err(HydrogenBondError::InvalidOptions)
    }
}

/// The minimum-image box for a periodic request, refusing absent or placeholder
/// cells.
pub(crate) fn periodic_box(
    structure: &Structure,
    options: HydrogenBondOptions,
    policy: HydrogenBondPolicy,
) -> Result<Option<PeriodicBox>, HydrogenBondError> {
    if !options.periodic {
        return Ok(None);
    }
    let cell = structure
        .data()
        .cell
        .ok_or(HydrogenBondError::MissingCell)?;
    if cell.is_placeholder() && !policy.allow_placeholder_cell {
        return Err(HydrogenBondError::MissingCell);
    }
    Ok(Some(PeriodicBox::from_cell(cell)?))
}

/// Whether two atoms can coexist in one conformation: a blank alternate
/// location fits any other, and two labels must be equal.
pub(crate) fn altlocs_compatible(structure: &Structure, first: u32, second: u32) -> bool {
    let label = |atom: u32| {
        structure
            .data()
            .atom(AtomIndex::new(atom))
            .and_then(AtomRef::alt_label)
    };
    match (label(first), label(second)) {
        (Some(a), Some(b)) => a == b,
        _ => true,
    }
}

/// Donor–acceptor heavy-atom pairs within the distance cutoff, with no
/// hydrogen geometry. Pairs are oriented `(donor, acceptor)`, sorted and
/// deduplicated; bonded and altloc-incompatible pairs are removed.
pub(crate) fn heavy_atom_pairs(
    structure: &Structure,
    options: HydrogenBondOptions,
    policy: HydrogenBondPolicy,
    context: &ExecutionContext,
) -> Result<Vec<(u32, u32)>, HydrogenBondError> {
    validate_options(options, policy)?;
    let donors = role_selection(structure, molframe_core::HBOND_DONOR_ANNOTATION)?;
    let acceptors = role_selection(structure, molframe_core::HBOND_ACCEPTOR_ANNOTATION)?;
    let periodic_box = periodic_box(structure, options, policy)?;
    let bonds = &structure.data().bonds;
    let adjacency = bonds
        .is_available()
        .then(|| bonds.adjacency(structure.atom_count()));
    let query = PairQuery {
        positions: structure.positions(),
        left: &donors,
        right: &acceptors,
        cutoff: options.maximum_donor_acceptor_distance,
        options: SpatialSearchOptions::with_backend(options.backend),
        periodic: periodic_box.as_ref(),
        context,
    };
    let parts =
        reduce_pairs_within_unsorted(&query, Vec::new, |out: &mut Vec<(u32, u32)>, pair| {
            for (donor, acceptor) in orientations(pair.first, pair.second, &donors, &acceptors) {
                let bonded = adjacency.is_some_and(|adjacency| {
                    adjacency
                        .neighbours(AtomIndex::new(donor))
                        .binary_search(&AtomIndex::new(acceptor))
                        .is_ok()
                });
                if donor != acceptor && !bonded && altlocs_compatible(structure, donor, acceptor) {
                    out.push((donor, acceptor));
                }
            }
        })?;
    let mut pairs: Vec<(u32, u32)> = parts.into_iter().flatten().collect();
    pairs.sort_unstable();
    pairs.dedup();
    Ok(pairs)
}

fn role_selection(structure: &Structure, name: &str) -> Result<AtomSelection, HydrogenBondError> {
    let Some(AtomAnnotation::Boolean(column)) = structure.annotations().get(name) else {
        return Err(HydrogenBondError::MissingChemistry);
    };
    if column.len() != structure.atom_count() {
        return Err(HydrogenBondError::MissingChemistry);
    }
    Ok(AtomSelection::from_sorted(
        (0..structure.atom_count())
            .filter(|atom| matches!(column.get(*atom), Some((true, Presence::Present))))
            .collect(),
    ))
}

fn orientations(
    first: u32,
    second: u32,
    donors: &AtomSelection,
    acceptors: &AtomSelection,
) -> impl Iterator<Item = (u32, u32)> {
    [
        (donors.contains(first) && acceptors.contains(second)).then_some((first, second)),
        (donors.contains(second) && acceptors.contains(first)).then_some((second, first)),
    ]
    .into_iter()
    .flatten()
}

fn donor_hydrogens<'a>(
    structure: &'a Structure,
    adjacency: &'a molframe_core::BondAdjacency,
    donor: u32,
) -> impl Iterator<Item = u32> + 'a {
    adjacency
        .neighbours(AtomIndex::new(donor))
        .iter()
        .filter_map(|atom| {
            structure
                .data()
                .atom(*atom)
                .and_then(AtomRef::element)
                .is_some_and(molframe_core::Element::is_hydrogen)
                .then_some(atom.get())
        })
}

fn measure(
    structure: &Structure,
    donor: u32,
    hydrogen: u32,
    acceptor: u32,
    periodic: Option<&PeriodicBox>,
) -> Option<HydrogenBond> {
    let positions = structure.positions();
    let donor_position = *positions.get(donor as usize)?;
    let hydrogen_position = *positions.get(hydrogen as usize)?;
    let acceptor_position = *positions.get(acceptor as usize)?;
    let (hydrogen_image, acceptor_image) = match periodic {
        Some(periodic) => (
            translated(
                donor_position,
                periodic.displacement(donor_position, hydrogen_position),
            ),
            translated(
                donor_position,
                periodic.displacement(donor_position, acceptor_position),
            ),
        ),
        None => (hydrogen_position, acceptor_position),
    };
    Some(HydrogenBond {
        donor: AtomIndex::new(donor),
        hydrogen: AtomIndex::new(hydrogen),
        acceptor: AtomIndex::new(acceptor),
        donor_acceptor_distance: f64_to_f32(molframe_geom::distance(
            donor_position,
            acceptor_image,
        )),
        hydrogen_acceptor_distance: f64_to_f32(molframe_geom::distance(
            hydrogen_image,
            acceptor_image,
        )),
        angle_degrees: molframe_geom::angle(donor_position, hydrogen_image, acceptor_image)?
            .to_degrees(),
    })
}

fn translated(origin: [f32; 3], displacement: [f32; 3]) -> [f32; 3] {
    std::array::from_fn(|axis| origin[axis] + displacement[axis])
}

#[cfg(test)]
#[path = "hbond_tests.rs"]
mod tests;
