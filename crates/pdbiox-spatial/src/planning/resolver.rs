//! Query-language spatial execution over a structure snapshot.

use crate::brute::finite;
use crate::{
    CellList, KdTree, NeighborPair, PeriodicBox, SpatialBackend, SpatialSearchOptions,
    pairs_within_with_options,
};
use pdbiox_core::contract::AnalysisPolicy;
use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::selection::AtomSelection;
use pdbiox_core::structure::Structure;
use pdbiox_query::{GeometricRequest, SpatialResolver};
use std::collections::HashMap;
use std::sync::{Arc, Mutex, MutexGuard};

#[path = "resolver_helpers.rs"]
mod helpers;
use helpers::{
    atom_out_of_bounds, backend_key, cacheable_backend, collect_within, index_cutoff_key,
    mark_pair_matches, periodic_box, radial_squared, resolve_beyond_request,
    resolve_within_request, select_iso_layer, spatial_diagnostic, update_nearest_distances,
    validate_shell, z_matches,
};

/// A spatial resolver bound to one immutable structure snapshot.
#[derive(Debug)]
pub struct StructureSpatial<'a> {
    structure: &'a Structure,
    options: SpatialSearchOptions,
    periodic: Option<PeriodicBox>,
    cache: Mutex<HashMap<IndexKey, Arc<CachedIndex<'a>>>>,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
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
        Self::new_with_options(
            structure,
            policy,
            SpatialSearchOptions::with_backend(backend),
        )
    }

    /// Creates a resolver under explicit spatial planning options.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic for invalid options or periodic cell data.
    pub fn new_with_options(
        structure: &'a Structure,
        policy: &AnalysisPolicy,
        options: SpatialSearchOptions,
    ) -> Result<Self, Diagnostic> {
        options.validate().map_err(spatial_diagnostic)?;
        let periodic = periodic_box(structure, policy)?;

        Ok(Self {
            structure,
            options,
            periodic,
            cache: Mutex::new(HashMap::new()),
        })
    }

    /// Number of structure-bound indices currently retained for reuse.
    #[must_use]
    pub fn cached_index_count(&self) -> usize {
        self.lock_cache().len()
    }

    /// Returns this structure snapshot's coordinate slice.
    ///
    /// Runtime and auxiliary space are `O(1)`.
    #[inline]
    fn positions(&self) -> &[[f32; 3]] {
        self.structure.positions()
    }

    /// Resolves one atom position or returns the query out-of-bounds diagnostic.
    fn position(&self, atom: u32) -> Result<[f32; 3], Diagnostic> {
        let index = usize::try_from(atom).map_err(|_| atom_out_of_bounds(atom))?;
        self.positions()
            .get(index)
            .copied()
            .ok_or_else(|| atom_out_of_bounds(atom))
    }

    /// Selects atoms inside a radial shell around `centre`.
    ///
    /// When `z` is provided, radial distance is measured in XY and the
    /// displacement must also lie inside the requested Z interval.
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
        let cylindrical = z.is_some();
        let mut selected = Vec::new();

        for atom in universe {
            let position = self.position(atom)?;

            if !finite(position) {
                continue;
            }

            let displacement = self.displacement(centre, position);

            if !z_matches(displacement[2], z) {
                continue;
            }

            let squared = radial_squared(displacement, cylindrical);

            if squared >= inner_squared && squared <= outer_squared {
                selected.push(atom);
            }
        }

        Ok(AtomSelection::from_sorted(selected))
    }

    /// Computes periodic or ordinary displacement from `left` to `right`.
    ///
    /// Runtime and auxiliary space are `O(1)`.
    #[inline]
    fn displacement(&self, left: [f32; 3], right: [f32; 3]) -> [f32; 3] {
        match &self.periodic {
            Some(periodic) => periodic.displacement(left, right),
            None => [right[0] - left[0], right[1] - left[1], right[2] - left[2]],
        }
    }

    /// Computes the arithmetic centre of finite atoms in `target`.
    ///
    /// # Errors
    ///
    /// Returns `E5003` for an empty/no-finite target and `E6009` for invalid
    /// atom indices.
    fn centre(&self, target: &AtomSelection) -> Result<[f32; 3], Diagnostic> {
        if target.is_empty() {
            return Err(Diagnostic::new(Code::E5003));
        }

        let mut sum = [0.0f64; 3];
        let mut count = 0u32;

        for atom in target {
            let position = self.position(atom)?;

            if !finite(position) {
                continue;
            }

            sum[0] += f64::from(position[0]);
            sum[1] += f64::from(position[1]);
            sum[2] += f64::from(position[2]);
            count += 1;
        }

        if count == 0 {
            return Err(Diagnostic::new(Code::E5003));
        }

        let divisor = f64::from(count);

        let convert = |value| {
            crate::numeric::finite_f32(value / divisor).ok_or_else(|| Diagnostic::new(Code::E5003))
        };
        Ok([convert(sum[0])?, convert(sum[1])?, convert(sum[2])?])
    }

    /// Resolves a shell around the centre of `target`.
    ///
    /// This helper centralises sphere/cylinder centre calculation.
    fn radial_from_target(
        &self,
        universe: &AtomSelection,
        target: &AtomSelection,
        inner: f32,
        outer: f32,
        z: Option<(f32, f32)>,
    ) -> Result<AtomSelection, Diagnostic> {
        let centre = self.centre(target)?;
        self.radial(universe, centre, inner, outer, z)
    }

    /// Selects atoms whose nearest target distance lies inside `[inner, outer]`.
    ///
    /// A dense nearest-distance buffer gives `O(P + U)` processing after the
    /// spatial pair search and avoids tree-map operations per pair.
    fn iso_layer(
        &self,
        universe: &AtomSelection,
        target: &AtomSelection,
        inner: f32,
        outer: f32,
    ) -> Result<AtomSelection, Diagnostic> {
        validate_shell(inner, outer)?;

        let pairs = self.pairs(universe, target, outer)?;
        let mut nearest = vec![f32::INFINITY; self.positions().len()];

        update_nearest_distances(&mut nearest, universe, target, pairs)?;

        select_iso_layer(universe, target, &nearest, inner * inner, outer * outer)
    }

    /// Selects query atoms within `cutoff` of at least one target.
    ///
    /// A dense bit vector removes the previous sort/dedup pass.
    fn within(
        &self,
        query: &AtomSelection,
        target: &AtomSelection,
        cutoff: f32,
    ) -> Result<AtomSelection, Diagnostic> {
        let pairs = self.pairs(query, target, cutoff)?;
        let mut matched = vec![false; self.positions().len()];

        mark_pair_matches(&mut matched, query, target, pairs)?;

        collect_within(query, target, &matched)
    }

    /// Dispatches pair search through periodic, cached or uncached backends.
    ///
    /// Non-periodic cell/k-d indices are cached for repeated geometric
    /// predicates over the same target selection.
    fn pairs(
        &self,
        query: &AtomSelection,
        target: &AtomSelection,
        cutoff: f32,
    ) -> Result<Vec<NeighborPair>, Diagnostic> {
        if self.periodic.is_some() {
            return pairs_within_with_options(
                self.positions(),
                query,
                target,
                cutoff,
                self.options,
                self.periodic.as_ref(),
            )
            .map_err(spatial_diagnostic);
        }

        let backend = cacheable_backend(self.options, query.len(), target.len())
            .map_err(spatial_diagnostic)?;

        if matches!(backend, SpatialBackend::CellList | SpatialBackend::KdTree) {
            self.cached_pairs(query, target, cutoff, backend)
        } else {
            let options = SpatialSearchOptions {
                backend,
                ..self.options
            };
            pairs_within_with_options(self.positions(), query, target, cutoff, options, None)
                .map_err(spatial_diagnostic)
        }
    }

    /// Executes a query using a cached `CellList` or `KdTree`.
    ///
    /// The returned `Arc` allows the mutex to be released before the expensive
    /// spatial query itself is executed.
    fn cached_pairs(
        &self,
        query: &AtomSelection,
        target: &AtomSelection,
        cutoff: f32,
        backend: SpatialBackend,
    ) -> Result<Vec<NeighborPair>, Diagnostic> {
        let selection: Vec<u32> = target.iter().collect();

        let key = IndexKey {
            generation: self.structure.generation().get(),
            selection: selection.into_boxed_slice(),
            cutoff: index_cutoff_key(backend, cutoff),
            backend: backend_key(backend),
        };

        let index = self.cached_index(key, backend, cutoff)?;

        let query: Vec<u32> = query.iter().collect();

        match index.as_ref() {
            CachedIndex::Cell(index) => index.pairs(&query, cutoff).map_err(spatial_diagnostic),
            CachedIndex::Kd(index) => index.pairs(&query, cutoff).map_err(spatial_diagnostic),
        }
    }

    /// Returns or constructs one cached structure-bound spatial index.
    ///
    /// Cache lookup is expected `O(1)`. An `Arc` is cloned so the mutex is not
    /// held while performing the actual query.
    fn cached_index(
        &self,
        key: IndexKey,
        backend: SpatialBackend,
        cutoff: f32,
    ) -> Result<Arc<CachedIndex<'a>>, Diagnostic> {
        let mut cache = self.lock_cache();

        if let Some(index) = cache.get(&key) {
            return Ok(Arc::clone(index));
        }

        let index = Arc::new(self.build_cached_index(backend, &key.selection, cutoff)?);

        cache.insert(key, Arc::clone(&index));
        Ok(index)
    }

    /// Constructs one cacheable index over `targets`.
    ///
    /// K-d trees are cutoff-independent; `CellLists` retain their build cutoff.
    fn build_cached_index(
        &self,
        backend: SpatialBackend,
        targets: &[u32],
        cutoff: f32,
    ) -> Result<CachedIndex<'a>, Diagnostic> {
        let positions: &'a [[f32; 3]] = self.structure.positions();

        match backend {
            SpatialBackend::CellList => CellList::build_with_options(
                positions,
                targets,
                cutoff,
                None,
                self.options.cell_grid,
            )
            .map(CachedIndex::Cell)
            .map_err(spatial_diagnostic),
            SpatialBackend::KdTree => {
                KdTree::build_with_options(positions, targets, None, self.options.kd_periodic)
                    .map(CachedIndex::Kd)
                    .map_err(spatial_diagnostic)
            }
            _ => Err(Diagnostic::new(Code::E9001)),
        }
    }

    /// Acquires the index cache, recovering data from a poisoned mutex.
    ///
    /// Poison recovery preserves the original behavior and never panics.
    fn lock_cache(&self) -> MutexGuard<'_, HashMap<IndexKey, Arc<CachedIndex<'a>>>> {
        match self.cache.lock() {
            Ok(cache) => cache,
            Err(poisoned) => poisoned.into_inner(),
        }
    }
}

impl SpatialResolver for StructureSpatial<'_> {
    /// Resolves one query-language geometric request.
    ///
    /// Each variant delegates to a focused helper to keep dispatch complexity
    /// independent from the spatial algorithms themselves.
    fn resolve(&self, request: GeometricRequest<'_>) -> Result<AtomSelection, Diagnostic> {
        match request {
            GeometricRequest::Within {
                universe,
                target,
                radius,
                include_target,
            } => resolve_within_request(self, universe, target, radius, include_target),
            GeometricRequest::Beyond {
                universe,
                target,
                radius,
            } => resolve_beyond_request(self, universe, target, radius),
            GeometricRequest::SphereZone {
                universe,
                target,
                radius,
            } => self.radial_from_target(universe, target, 0.0, radius, None),
            GeometricRequest::SphereLayer {
                universe,
                target,
                inner,
                outer,
            } => self.radial_from_target(universe, target, inner, outer, None),
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
            } => self.radial_from_target(universe, target, 0.0, radius, Some((z_min, z_max))),
            GeometricRequest::CylinderLayer {
                universe,
                target,
                inner,
                outer,
                z_max,
                z_min,
            } => self.radial_from_target(universe, target, inner, outer, Some((z_min, z_max))),
            GeometricRequest::Point {
                universe,
                point,
                radius,
            } => self.radial(universe, point, 0.0, radius, None),
            _ => Err(Diagnostic::new(Code::E9001)),
        }
    }
}

#[cfg(test)]
#[path = "resolver_tests.rs"]
mod tests;
