//! Symmetry-aware crystal neighbour enumeration.

use crate::SymmetrySet;
use crate::crystal_images::relevant_images;
use pdbiox_core::selection::AtomSelection;
use pdbiox_core::{AtomIndex, Code, Diagnostic, ModelIndex, Structure};
use pdbiox_spatial::{SpatialBackend, pairs_within};
use std::collections::BTreeSet;

/// Default ceiling on candidate atom images examined by a crystal search.
pub const DEFAULT_CRYSTAL_IMAGE_LIMIT: usize = 10_000_000;

/// One unique contact from an asymmetric-unit atom to a crystal image.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CrystalNeighbor {
    /// Atom anchored in the deposited asymmetric unit.
    pub source_atom: AtomIndex,
    /// Source atom transformed into the neighbouring image.
    pub image_atom: AtomIndex,
    /// Position of the explicit symmetry representative.
    pub operation: usize,
    /// Primitive lattice translation after applying the representative.
    pub lattice: [i32; 3],
    /// Squared Cartesian distance in ångström².
    pub distance_squared: f64,
}

/// Finds unique atom pairs across explicit crystallographic images.
///
/// # Errors
///
/// Returns a diagnostic for missing crystal data, invalid cutoff/model, or
/// work above the default candidate-image limit.
pub fn crystal_neighbors(
    structure: &Structure,
    symmetry: &SymmetrySet,
    model: ModelIndex,
    cutoff: f64,
) -> Result<Vec<CrystalNeighbor>, Diagnostic> {
    crystal_neighbors_with_backend(
        structure,
        symmetry,
        model,
        cutoff,
        SpatialBackend::Auto,
        DEFAULT_CRYSTAL_IMAGE_LIMIT,
    )
}

/// Crystal neighbour search under an explicit candidate-image ceiling.
///
/// # Errors
///
/// Returns a registered diagnostic before exceeding `limit` candidates.
pub fn crystal_neighbors_with_limit(
    structure: &Structure,
    symmetry: &SymmetrySet,
    model: ModelIndex,
    cutoff: f64,
    limit: usize,
) -> Result<Vec<CrystalNeighbor>, Diagnostic> {
    crystal_neighbors_with_backend(
        structure,
        symmetry,
        model,
        cutoff,
        SpatialBackend::Auto,
        limit,
    )
}

/// Crystal neighbour search using a selected spatial backend.
///
/// # Errors
///
/// Returns the same diagnostics as [`crystal_neighbors_with_limit`].
pub fn crystal_neighbors_with_backend(
    structure: &Structure,
    symmetry: &SymmetrySet,
    model: ModelIndex,
    cutoff: f64,
    backend: SpatialBackend,
    limit: usize,
) -> Result<Vec<CrystalNeighbor>, Diagnostic> {
    let (sources, images) = relevant_images(structure, symmetry, model, cutoff, limit)?;
    let source_count = u32::try_from(sources.len()).map_err(|_| search_limit())?;
    let total_count = sources
        .len()
        .checked_add(images.len())
        .ok_or_else(search_limit)?;
    let total_count = u32::try_from(total_count).map_err(|_| search_limit())?;
    let mut positions = Vec::with_capacity(total_count as usize);
    positions.extend(sources.iter().map(|source| source.cartesian));
    positions.extend(images.iter().map(|image| image.cartesian));
    let left = AtomSelection::from_sorted((0..source_count).collect());
    let right = AtomSelection::from_sorted((source_count..total_count).collect());
    let pairs = pairs_within(&positions, &left, &right, cutoff as f32, backend, None)
        .map_err(pdbiox_spatial::SpatialError::into_diagnostic)?;
    let mut seen = BTreeSet::new();
    let mut neighbors = Vec::new();
    for pair in pairs {
        let source = sources.get(pair.first as usize).ok_or_else(invariant)?;
        let image_index = pair
            .second
            .checked_sub(source_count)
            .ok_or_else(invariant)?;
        let image = images.get(image_index as usize).ok_or_else(invariant)?;
        let distance_squared = f64::from(pair.distance_squared);
        if source.atom == image.atom && distance_squared <= f64::EPSILON {
            continue;
        }
        let forward = PairKey::new(source.atom, image.atom, image.operation, image.lattice)?;
        let (inverse_operation, inverse_lattice) =
            symmetry.inverse_image(image.operation, image.lattice)?;
        let reverse = PairKey::new(image.atom, source.atom, inverse_operation, inverse_lattice)?;
        let canonical = forward.min(reverse);
        if seen.insert(canonical) {
            neighbors.push(CrystalNeighbor {
                source_atom: AtomIndex::new(canonical.source),
                image_atom: AtomIndex::new(canonical.image),
                operation: canonical.operation as usize,
                lattice: canonical.lattice,
                distance_squared,
            });
        }
    }
    neighbors.sort_by(|left, right| {
        left.source_atom
            .cmp(&right.source_atom)
            .then(left.image_atom.cmp(&right.image_atom))
            .then(left.operation.cmp(&right.operation))
            .then(left.lattice.cmp(&right.lattice))
    });
    Ok(neighbors)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct PairKey {
    source: u32,
    image: u32,
    operation: u32,
    lattice: [i32; 3],
}

impl PairKey {
    fn new(
        source: AtomIndex,
        image: AtomIndex,
        operation: usize,
        lattice: [i32; 3],
    ) -> Result<Self, Diagnostic> {
        Ok(Self {
            source: source.get(),
            image: image.get(),
            operation: u32::try_from(operation).map_err(|_| search_limit())?,
            lattice,
        })
    }
}

fn search_limit() -> Diagnostic {
    Diagnostic::new(Code::E6017)
}

fn invariant() -> Diagnostic {
    Diagnostic::new(Code::E9001)
}

#[cfg(test)]
#[path = "crystal_tests.rs"]
mod tests;
