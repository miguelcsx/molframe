//! CCD-perceived hydrogen bonds with explicit hydrogen geometry.

use pdbiox_core::column::Presence;
use pdbiox_core::index::AtomIndex;
use pdbiox_core::selection::AtomSelection;
use pdbiox_core::structure::{AtomRef, Structure};
use pdbiox_core::{ExecutionContext, annotation::AtomAnnotation};
use pdbiox_spatial::{
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
    /// Periodic geometry was requested without a unit cell.
    #[error("periodic hydrogen bonds require a unit cell")]
    MissingCell,
    /// Spatial indexing rejected the request.
    #[error(transparent)]
    Spatial(#[from] SpatialError),
}

/// Finds hydrogen bonds from CCD roles, explicit hydrogens and angular geometry.
///
/// No heavy-atom-only direction or residue-name fallback is performed. Results
/// are sorted by `(donor, hydrogen, acceptor)`.
///
/// # Errors
///
/// Returns [`HydrogenBondError`] for invalid options, absent chemistry,
/// unavailable periodic cell data or a spatial-search failure.
pub fn hydrogen_bonds(
    structure: &Structure,
    options: HydrogenBondOptions,
    context: &ExecutionContext,
) -> Result<Vec<HydrogenBond>, HydrogenBondError> {
    validate_options(options)?;
    if !structure.data().bonds.is_available() {
        return Err(HydrogenBondError::MissingChemistry);
    }
    let donors = role_selection(structure, pdbiox_core::HBOND_DONOR_ANNOTATION)?;
    let acceptors = role_selection(structure, pdbiox_core::HBOND_ACCEPTOR_ANNOTATION)?;
    let cell = if options.periodic {
        Some(
            structure
                .data()
                .cell
                .ok_or(HydrogenBondError::MissingCell)?,
        )
    } else {
        None
    };
    let periodic_box = cell.map(PeriodicBox::from_cell).transpose()?;
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
    Ok(output)
}

fn validate_options(options: HydrogenBondOptions) -> Result<(), HydrogenBondError> {
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
    adjacency: &'a pdbiox_core::BondAdjacency,
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
                .is_some_and(pdbiox_core::Element::is_hydrogen)
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
        donor_acceptor_distance: f64_to_f32(pdbiox_geom::distance(donor_position, acceptor_image)),
        hydrogen_acceptor_distance: f64_to_f32(pdbiox_geom::distance(
            hydrogen_image,
            acceptor_image,
        )),
        angle_degrees: pdbiox_geom::angle(donor_position, hydrogen_image, acceptor_image)?
            .to_degrees(),
    })
}

fn translated(origin: [f32; 3], displacement: [f32; 3]) -> [f32; 3] {
    std::array::from_fn(|axis| origin[axis] + displacement[axis])
}

#[cfg(test)]
#[path = "hbond_tests.rs"]
mod tests;
