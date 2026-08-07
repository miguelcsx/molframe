//! Query-language spatial execution over a structure snapshot.

use crate::{CellList, KdTree, PeriodicBox, SpatialBackend, SpatialError, pairs_within};
use pdbiox_core::contract::{AnalysisPolicy, PeriodicPolicy};
use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::selection::AtomSelection;
use pdbiox_core::structure::Structure;
use pdbiox_query::{GeometricRequest, SpatialResolver};
use std::collections::BTreeMap;
use std::sync::Mutex;

/// A spatial resolver bound to one immutable structure snapshot.
#[derive(Debug)]
pub struct StructureSpatial<'a> {
    structure: &'a Structure,
    backend: SpatialBackend,
    periodic: Option<PeriodicBox>,
    cache: Mutex<BTreeMap<IndexKey, CachedIndex<'a>>>,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct IndexKey {
    generation: u64,
    selection: Box<[u32]>,
    cutoff: u32,
    backend: u8,
}

#[derive(Debug)]
enum CachedIndex<'a> {
    Cell(CellList<'a>),
    Kd(KdTree<'a>),
}

impl<'a> StructureSpatial<'a> {
    /// Creates a resolver under the policy's periodicity decision.
    ///
    /// # Errors
    ///
    /// Returns `E5004` when periodicity is requested without a valid unit cell.
    pub fn new(
        structure: &'a Structure,
        policy: &AnalysisPolicy,
        backend: SpatialBackend,
    ) -> Result<Self, Diagnostic> {
        let periodic = match policy.periodic {
            PeriodicPolicy::None => None,
            PeriodicPolicy::Pbc | PeriodicPolicy::MinimumImage => {
                let Some(cell) = structure.data().cell else {
                    return Err(Diagnostic::new(Code::E5004));
                };
                Some(PeriodicBox::from_cell(cell).map_err(spatial_diagnostic)?)
            }
            _ => return Err(Diagnostic::new(Code::E6004)),
        };
        Ok(Self {
            structure,
            backend,
            periodic,
            cache: Mutex::new(BTreeMap::new()),
        })
    }

    /// Number of structure-bound indices currently retained for reuse.
    #[must_use]
    pub fn cached_index_count(&self) -> usize {
        match self.cache.lock() {
            Ok(cache) => cache.len(),
            Err(poisoned) => poisoned.into_inner().len(),
        }
    }

    fn positions(&self) -> &[[f32; 3]] {
        self.structure.positions()
    }

    fn radial(
        &self,
        universe: &AtomSelection,
        centre: [f32; 3],
        inner: f32,
        outer: f32,
        z: Option<(f32, f32)>,
    ) -> Result<AtomSelection, Diagnostic> {
        validate_shell(inner, outer)?;
        let inner_squared = inner * inner;
        let outer_squared = outer * outer;
        let mut selected = Vec::new();
        for atom in universe {
            let Some(position) = self.positions().get(atom as usize).copied() else {
                return Err(Diagnostic::new(Code::E6009).with_context("atom", atom.to_string()));
            };
            if !position.iter().all(|component| component.is_finite()) {
                continue;
            }
            let displacement = match &self.periodic {
                Some(periodic) => periodic.displacement(centre, position),
                None => [
                    position[0] - centre[0],
                    position[1] - centre[1],
                    position[2] - centre[2],
                ],
            };
            if let Some((z_min, z_max)) = z
                && (displacement[2] < z_min || displacement[2] > z_max)
            {
                continue;
            }
            let axes = if z.is_some() { 2 } else { 3 };
            let squared: f32 = displacement[..axes]
                .iter()
                .map(|component| component * component)
                .sum();
            if squared >= inner_squared && squared <= outer_squared {
                selected.push(atom);
            }
        }
        Ok(AtomSelection::from_sorted(selected))
    }

    fn centre(&self, target: &AtomSelection) -> Result<[f32; 3], Diagnostic> {
        if target.is_empty() {
            return Err(Diagnostic::new(Code::E5003));
        }
        let mut sum = [0.0f64; 3];
        let mut count = 0u32;
        for atom in target {
            let Some(position) = self.positions().get(atom as usize) else {
                return Err(Diagnostic::new(Code::E6009).with_context("atom", atom.to_string()));
            };
            if position.iter().all(|component| component.is_finite()) {
                for axis in 0..3 {
                    sum[axis] += f64::from(position[axis]);
                }
                count += 1;
            }
        }
        if count == 0 {
            return Err(Diagnostic::new(Code::E5003));
        }
        Ok(sum.map(|value| (value / f64::from(count)) as f32))
    }

    fn iso_layer(
        &self,
        universe: &AtomSelection,
        target: &AtomSelection,
        inner: f32,
        outer: f32,
    ) -> Result<AtomSelection, Diagnostic> {
        validate_shell(inner, outer)?;
        let pairs = self.pairs(universe, target, outer)?;
        let mut nearest: BTreeMap<u32, f32> = BTreeMap::new();
        let overlap = universe.intersect(target);
        for atom in &overlap {
            nearest.insert(atom, 0.0);
        }
        for pair in pairs {
            if universe.contains(pair.first) && target.contains(pair.second) {
                keep_nearest(&mut nearest, pair.first, pair.distance_squared);
            }
            if universe.contains(pair.second) && target.contains(pair.first) {
                keep_nearest(&mut nearest, pair.second, pair.distance_squared);
            }
        }
        let inner_squared = inner * inner;
        Ok(AtomSelection::from_sorted(
            nearest
                .into_iter()
                .filter_map(|(atom, squared)| (squared >= inner_squared).then_some(atom))
                .collect(),
        ))
    }

    fn within(
        &self,
        query: &AtomSelection,
        target: &AtomSelection,
        cutoff: f32,
    ) -> Result<AtomSelection, Diagnostic> {
        let pairs = self.pairs(query, target, cutoff)?;
        let mut selected: Vec<u32> = query.intersect(target).iter().collect();
        for pair in pairs {
            if query.contains(pair.first) && target.contains(pair.second) {
                selected.push(pair.first);
            }
            if query.contains(pair.second) && target.contains(pair.first) {
                selected.push(pair.second);
            }
        }
        selected.sort_unstable();
        selected.dedup();
        Ok(AtomSelection::from_sorted(selected))
    }

    fn pairs(
        &self,
        query: &AtomSelection,
        target: &AtomSelection,
        cutoff: f32,
    ) -> Result<Vec<crate::NeighborPair>, Diagnostic> {
        if self.periodic.is_some() {
            return pairs_within(
                self.positions(),
                query,
                target,
                cutoff,
                self.backend,
                self.periodic.as_ref(),
            )
            .map_err(spatial_diagnostic);
        }
        let backend = cacheable_backend(self.backend, query.len(), target.len());
        if !matches!(backend, SpatialBackend::CellList | SpatialBackend::KdTree) {
            return pairs_within(self.positions(), query, target, cutoff, backend, None)
                .map_err(spatial_diagnostic);
        }
        let targets: Vec<u32> = target.iter().collect();
        let key = IndexKey {
            generation: self.structure.generation().get(),
            selection: targets.clone().into_boxed_slice(),
            cutoff: cutoff.to_bits(),
            backend: backend_key(backend),
        };
        let mut cache = match self.cache.lock() {
            Ok(cache) => cache,
            Err(poisoned) => poisoned.into_inner(),
        };
        if !cache.contains_key(&key) {
            let structure: &'a Structure = self.structure;
            let index = match backend {
                SpatialBackend::CellList => CachedIndex::Cell(
                    CellList::build(structure.positions(), &targets, cutoff, None)
                        .map_err(spatial_diagnostic)?,
                ),
                SpatialBackend::KdTree => CachedIndex::Kd(
                    KdTree::build(structure.positions(), &targets, None)
                        .map_err(spatial_diagnostic)?,
                ),
                _ => return Err(Diagnostic::new(Code::E9001)),
            };
            cache.insert(key.clone(), index);
        }
        let query: Vec<u32> = query.iter().collect();
        match cache.get(&key) {
            Some(CachedIndex::Cell(index)) => {
                index.pairs(&query, cutoff).map_err(spatial_diagnostic)
            }
            Some(CachedIndex::Kd(index)) => index.pairs(&query, cutoff).map_err(spatial_diagnostic),
            None => Err(Diagnostic::new(Code::E9001)),
        }
    }
}

impl SpatialResolver for StructureSpatial<'_> {
    fn resolve(&self, request: GeometricRequest<'_>) -> Result<AtomSelection, Diagnostic> {
        match request {
            GeometricRequest::Within {
                universe,
                target,
                radius,
                include_target,
            } => {
                let selected = self.within(universe, target, radius)?;
                if include_target {
                    Ok(selected)
                } else {
                    Ok(selected.difference(target))
                }
            }
            GeometricRequest::Beyond {
                universe,
                target,
                radius,
            } => {
                let close = self.within(universe, target, radius)?;
                Ok(universe.difference(&close))
            }
            GeometricRequest::SphereZone {
                universe,
                target,
                radius,
            } => self.radial(universe, self.centre(target)?, 0.0, radius, None),
            GeometricRequest::SphereLayer {
                universe,
                target,
                inner,
                outer,
            } => self.radial(universe, self.centre(target)?, inner, outer, None),
            GeometricRequest::IsoLayer {
                universe,
                target,
                inner,
                outer,
            } => self.iso_layer(universe, target, inner, outer),
            GeometricRequest::CylinderZone {
                universe,
                target,
                radius,
                z_max,
                z_min,
            } => self.radial(
                universe,
                self.centre(target)?,
                0.0,
                radius,
                Some((z_min, z_max)),
            ),
            GeometricRequest::CylinderLayer {
                universe,
                target,
                inner,
                outer,
                z_max,
                z_min,
            } => self.radial(
                universe,
                self.centre(target)?,
                inner,
                outer,
                Some((z_min, z_max)),
            ),
            GeometricRequest::Point {
                universe,
                point,
                radius,
            } => self.radial(universe, point, 0.0, radius, None),
            _ => Err(Diagnostic::new(Code::E9001)),
        }
    }
}

fn cacheable_backend(backend: SpatialBackend, left: u32, right: u32) -> SpatialBackend {
    match backend {
        SpatialBackend::Auto if u64::from(left) * u64::from(right) > 250_000 => {
            if right >= 20_000 && left < right / 8 {
                SpatialBackend::KdTree
            } else {
                SpatialBackend::CellList
            }
        }
        SpatialBackend::Auto => SpatialBackend::BruteForce,
        other => other,
    }
}

fn backend_key(backend: SpatialBackend) -> u8 {
    match backend {
        SpatialBackend::CellList => 1,
        SpatialBackend::KdTree => 2,
        SpatialBackend::BruteForce | SpatialBackend::NeighborList | SpatialBackend::Auto => 0,
    }
}

fn validate_shell(inner: f32, outer: f32) -> Result<(), Diagnostic> {
    if inner.is_finite() && outer.is_finite() && inner >= 0.0 && outer >= inner {
        Ok(())
    } else {
        Err(Diagnostic::new(Code::E4002).with_context("geometry", "invalid radius interval"))
    }
}

fn keep_nearest(nearest: &mut BTreeMap<u32, f32>, atom: u32, squared: f32) {
    nearest
        .entry(atom)
        .and_modify(|current| *current = current.min(squared))
        .or_insert(squared);
}

fn spatial_diagnostic(error: SpatialError) -> Diagnostic {
    error.into_diagnostic()
}

#[cfg(test)]
#[path = "resolver_tests.rs"]
mod tests;
