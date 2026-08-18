//! Focused request, distance and cache helpers for structure-bound queries.

use super::StructureSpatial;
use crate::{NeighborPair, PeriodicBox, SpatialBackend, SpatialError, SpatialSearchOptions};
use pdbiox_core::contract::{AnalysisPolicy, PeriodicPolicy};
use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::selection::AtomSelection;
use pdbiox_core::structure::Structure;

pub(super) fn periodic_box(
    structure: &Structure,
    policy: &AnalysisPolicy,
) -> Result<Option<PeriodicBox>, Diagnostic> {
    match policy.periodic {
        PeriodicPolicy::None => Ok(None),
        PeriodicPolicy::Pbc | PeriodicPolicy::MinimumImage => {
            let Some(cell) = structure.data().cell else {
                return Err(Diagnostic::new(Code::E5004));
            };
            PeriodicBox::from_cell(cell)
                .map(Some)
                .map_err(spatial_diagnostic)
        }
        _ => Err(Diagnostic::new(Code::E6004)),
    }
}

pub(super) fn resolve_within_request(
    spatial: &StructureSpatial<'_>,
    universe: &AtomSelection,
    target: &AtomSelection,
    radius: f32,
    include_target: bool,
) -> Result<AtomSelection, Diagnostic> {
    let selected = spatial.within(universe, target, radius)?;
    if include_target {
        Ok(selected)
    } else {
        Ok(selected.difference(target))
    }
}

pub(super) fn resolve_beyond_request(
    spatial: &StructureSpatial<'_>,
    universe: &AtomSelection,
    target: &AtomSelection,
    radius: f32,
) -> Result<AtomSelection, Diagnostic> {
    let close = spatial.within(universe, target, radius)?;
    Ok(universe.difference(&close))
}

pub(super) fn update_nearest_distances(
    nearest: &mut [f32],
    universe: &AtomSelection,
    target: &AtomSelection,
    pairs: Vec<NeighborPair>,
) -> Result<(), Diagnostic> {
    for pair in pairs {
        if universe.contains(pair.first) && target.contains(pair.second) {
            keep_nearest(nearest, pair.first, pair.distance_squared)?;
        }
        if universe.contains(pair.second) && target.contains(pair.first) {
            keep_nearest(nearest, pair.second, pair.distance_squared)?;
        }
    }
    Ok(())
}

pub(super) fn select_iso_layer(
    universe: &AtomSelection,
    target: &AtomSelection,
    nearest: &[f32],
    inner_squared: f32,
    outer_squared: f32,
) -> Result<AtomSelection, Diagnostic> {
    let mut selected = Vec::new();
    for atom in universe {
        let index = usize::try_from(atom).map_err(|_| atom_out_of_bounds(atom))?;
        let squared = if target.contains(atom) {
            0.0
        } else {
            nearest
                .get(index)
                .copied()
                .ok_or_else(|| atom_out_of_bounds(atom))?
        };
        if squared >= inner_squared && squared <= outer_squared {
            selected.push(atom);
        }
    }
    Ok(AtomSelection::from_sorted(selected))
}

pub(super) fn mark_pair_matches(
    matched: &mut [bool],
    query: &AtomSelection,
    target: &AtomSelection,
    pairs: Vec<NeighborPair>,
) -> Result<(), Diagnostic> {
    for pair in pairs {
        if query.contains(pair.first) && target.contains(pair.second) {
            mark_atom(matched, pair.first)?;
        }
        if query.contains(pair.second) && target.contains(pair.first) {
            mark_atom(matched, pair.second)?;
        }
    }
    Ok(())
}

pub(super) fn collect_within(
    query: &AtomSelection,
    target: &AtomSelection,
    matched: &[bool],
) -> Result<AtomSelection, Diagnostic> {
    let mut selected = Vec::new();
    for atom in query {
        let index = usize::try_from(atom).map_err(|_| atom_out_of_bounds(atom))?;
        let spatial_match = matched
            .get(index)
            .copied()
            .ok_or_else(|| atom_out_of_bounds(atom))?;
        if spatial_match || target.contains(atom) {
            selected.push(atom);
        }
    }
    Ok(AtomSelection::from_sorted(selected))
}

pub(super) fn collect_selection_indices(selection: &AtomSelection) -> Vec<u32> {
    let Ok(capacity) = usize::try_from(selection.len()) else {
        return selection.iter().collect();
    };
    let mut indices = Vec::with_capacity(capacity);
    indices.extend(selection.iter());
    indices
}

fn mark_atom(matched: &mut [bool], atom: u32) -> Result<(), Diagnostic> {
    let index = usize::try_from(atom).map_err(|_| atom_out_of_bounds(atom))?;
    let Some(slot) = matched.get_mut(index) else {
        return Err(atom_out_of_bounds(atom));
    };
    *slot = true;
    Ok(())
}

fn keep_nearest(nearest: &mut [f32], atom: u32, squared: f32) -> Result<(), Diagnostic> {
    let index = usize::try_from(atom).map_err(|_| atom_out_of_bounds(atom))?;
    let Some(current) = nearest.get_mut(index) else {
        return Err(atom_out_of_bounds(atom));
    };
    *current = current.min(squared);
    Ok(())
}

#[inline]
pub(super) fn z_matches(displacement_z: f32, z: Option<(f32, f32)>) -> bool {
    match z {
        Some((z_min, z_max)) => displacement_z >= z_min && displacement_z <= z_max,
        None => true,
    }
}

#[inline]
pub(super) fn radial_squared(displacement: [f32; 3], cylindrical: bool) -> f32 {
    let xy = displacement[0] * displacement[0] + displacement[1] * displacement[1];
    if cylindrical {
        xy
    } else {
        xy + displacement[2] * displacement[2]
    }
}

pub(super) fn cacheable_backend(
    options: SpatialSearchOptions,
    left: u64,
    right: u64,
) -> Result<SpatialBackend, SpatialError> {
    let left = usize::try_from(left).map_err(|_| SpatialError::NumericRangeExceeded)?;
    let right = usize::try_from(right).map_err(|_| SpatialError::NumericRangeExceeded)?;
    options
        .plan(left, right, false, 0.0)
        .map(|plan| plan.backend)
}

#[inline]
pub(super) fn backend_key(backend: SpatialBackend) -> u8 {
    match backend {
        SpatialBackend::CellList => 1,
        SpatialBackend::KdTree => 2,
        SpatialBackend::BruteForce | SpatialBackend::NeighborList | SpatialBackend::Auto => 0,
    }
}

#[inline]
pub(super) fn index_cutoff_key(backend: SpatialBackend, cutoff: f32) -> u32 {
    match backend {
        SpatialBackend::CellList => cutoff.to_bits(),
        SpatialBackend::KdTree
        | SpatialBackend::BruteForce
        | SpatialBackend::NeighborList
        | SpatialBackend::Auto => 0,
    }
}

pub(super) fn validate_shell(inner: f32, outer: f32) -> Result<(), Diagnostic> {
    if inner.is_finite() && outer.is_finite() && inner >= 0.0 && outer >= inner {
        Ok(())
    } else {
        Err(Diagnostic::new(Code::E4002).with_context("geometry", "invalid radius interval"))
    }
}

pub(super) fn atom_out_of_bounds(atom: u32) -> Diagnostic {
    Diagnostic::new(Code::E6009).with_context("atom", atom.to_string())
}

pub(super) fn spatial_diagnostic(error: SpatialError) -> Diagnostic {
    error.into_diagnostic()
}
