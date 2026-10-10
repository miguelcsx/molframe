//! Dictionary-free canonical pairing geometry, O(residues + local candidates).

use crate::{CellList, SpatialError};
use molframe_chem::base_pair::{self as edges, CanonicalBase, WatsonCrickGeometry};
use molframe_core::{ResidueIndex, Structure};
use num_traits::ToPrimitive;

/// Canonical edge contacts and coplanar sugar geometry supporting a pair.
/// This is geometric evidence; it does not report perceived hydrogen bonds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GeometricBasePair {
    /// Lower-indexed nucleotide residue.
    pub first: ResidueIndex,
    /// Higher-indexed nucleotide residue.
    pub second: ResidueIndex,
    /// Number of canonical edge contacts within the requested distance.
    pub contact_count: usize,
    /// Shortest canonical edge contact, in ångström.
    pub closest_distance: f32,
}

/// Finds standard A–T/U and G–C pairs without a component dictionary.
///
/// Requires positioned C1', C2, C4 and C6 atoms, reference-plane agreement,
/// and at least two complementary edge contacts. Unknown or modified component
/// names are skipped; dictionary-backed chemistry is available through
/// the `molframe::analysis` domain. Positions are read from the first model.
///
/// # Errors
/// Returns a spatial error for non-finite or non-positive distance limits.
pub fn standard_base_pair_geometry(
    structure: &Structure,
    geometry: WatsonCrickGeometry,
    maximum_contact_distance: f32,
) -> Result<Vec<GeometricBasePair>, SpatialError> {
    if !maximum_contact_distance.is_finite()
        || maximum_contact_distance <= 0.0
        || !geometry.minimum_c1_distance.is_finite()
        || geometry.minimum_c1_distance < 0.0
        || !geometry.maximum_c1_distance.is_finite()
        || geometry.maximum_c1_distance <= geometry.minimum_c1_distance
        || !geometry.maximum_plane_angle_degrees.is_finite()
        || !(0.0..=90.0).contains(&geometry.maximum_plane_angle_degrees)
    {
        return Err(SpatialError::InvalidCutoff);
    }
    let mut residues = Vec::new();
    let mut positions = Vec::new();
    for residue in structure.data().residues() {
        let Some(base) = residue.name().and_then(standard_base) else {
            continue;
        };
        let Some(position) = residue
            .atom("C1'")
            .and_then(molframe_core::structure::AtomRef::position)
        else {
            continue;
        };
        residues.push((residue, base));
        positions.push(position);
    }
    let rows: Vec<_> = (0..positions.len())
        .map(u32::try_from)
        .collect::<Result<_, _>>()
        .map_err(|_| SpatialError::InvalidCutoff)?;
    let Some(cutoff) = geometry.maximum_c1_distance.to_f32() else {
        return Err(SpatialError::InvalidCutoff);
    };
    let index = CellList::build(&positions, &rows, cutoff, None)?;
    let mut pairs = Vec::new();
    for (row, &(left, left_base)) in residues.iter().enumerate() {
        index.for_each_neighbor(positions[row], cutoff, |other, _| {
            let other = other as usize;
            if other <= row {
                return;
            }
            let (right, right_base) = residues[other];
            if !left_base.complementary(right_base)
                || !edges::has_watson_crick_geometry(left, right, geometry)
            {
                return;
            }
            let mut count = 0;
            let mut closest = f32::INFINITY;
            for first in ["N1", "N2", "N3", "N4", "N6", "O2", "O4", "O6"] {
                for second in ["N1", "N2", "N3", "N4", "N6", "O2", "O4", "O6"] {
                    if !edges::watson_crick_contact(left_base, first, right_base, second) {
                        continue;
                    }
                    let (Some(a), Some(b)) = (
                        left.atom(first)
                            .and_then(molframe_core::structure::AtomRef::position),
                        right
                            .atom(second)
                            .and_then(molframe_core::structure::AtomRef::position),
                    ) else {
                        continue;
                    };
                    let squared: f32 = a.into_iter().zip(b).map(|(a, b)| (a - b).powi(2)).sum();
                    if squared <= maximum_contact_distance.powi(2) {
                        count += 1;
                        closest = closest.min(squared.sqrt());
                    }
                }
            }
            if count >= 2 {
                pairs.push(GeometricBasePair {
                    first: left.index(),
                    second: right.index(),
                    contact_count: count,
                    closest_distance: closest,
                });
            }
        })?;
    }
    pairs.sort_by_key(|pair| (pair.first, pair.second));
    Ok(pairs)
}

fn standard_base(name: &str) -> Option<CanonicalBase> {
    match name {
        "A" | "DA" => Some(CanonicalBase::Adenine),
        "C" | "DC" => Some(CanonicalBase::Cytosine),
        "G" | "DG" => Some(CanonicalBase::Guanine),
        "T" | "DT" => Some(CanonicalBase::Thymine),
        "U" | "DU" => Some(CanonicalBase::Uracil),
        _ => None,
    }
}

#[cfg(test)]
#[path = "base_pair_geometry_tests.rs"]
mod tests;
