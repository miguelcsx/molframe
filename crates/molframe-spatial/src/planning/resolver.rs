//! Query-language spatial execution over a structure snapshot.

use crate::brute::finite;
use crate::{
    CellList, KdTree, NeighborPair, PeriodicBox, SpatialBackend, SpatialSearchOptions,
    pairs_within_with_options,
};
use molframe_core::ExecutionContext;
use molframe_core::contract::AnalysisPolicy;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::selection::AtomSelection;
use molframe_core::structure::Structure;
use std::collections::{HashMap, hash_map::Entry};
use std::sync::{Arc, RwLock, RwLockReadGuard, RwLockWriteGuard};

#[path = "resolver_helpers.rs"]
mod helpers;
#[path = "resolver_dispatch.rs"]
mod resolver_dispatch;
#[path = "resolver_visit.rs"]
mod resolver_visit;
use helpers::{
    atom_out_of_bounds, backend_key, cacheable_backend, collect_selection_indices, collect_within,
    index_cutoff_key, mark_pair_match, periodic_box, radial_squared, select_iso_layer,
    spatial_diagnostic, update_nearest_distance, validate_shell, z_matches,
};

/// Maximum number of spatial indices retained by one resolver.
///
/// Query predicates can produce arbitrarily many distinct target selections
/// and cutoffs. Bounding retained indices prevents a long-lived resolver from
/// accumulating unbounded query-controlled cache state.
const MAX_CACHED_INDICES: usize = 64;

/// A spatial resolver bound to one immutable structure snapshot.
#[derive(Debug)]
pub struct StructureSpatial<'a> {
    structure: &'a Structure,
    context: &'a ExecutionContext,
    options: SpatialSearchOptions,
    periodic: Option<PeriodicBox>,
    cache: RwLock<HashMap<IndexKey, Arc<CachedIndex<'a>>>>,
}

/// Identifies a reusable spatial index within one immutable resolver.
///
/// Backend and cutoff precede the potentially large atom selection so equality
/// checks can reject incompatible plans before comparing selection contents.
#[derive(Debug, PartialEq, Eq, Hash)]
struct IndexKey {
    backend: u8,
    cutoff: u32,
    selection: Box<[u32]>,
}

/// A spatial index retained for reuse by compatible geometric predicates.
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
        context: &'a ExecutionContext,
    ) -> Result<Self, Diagnostic> {
        Self::new_with_options(
            structure,
            policy,
            SpatialSearchOptions::with_backend(backend),
            context,
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
        context: &'a ExecutionContext,
    ) -> Result<Self, Diagnostic> {
        options.validate().map_err(spatial_diagnostic)?;
        let periodic = periodic_box(structure, policy)?;

        Ok(Self {
            structure,
            context,
            options,
            periodic,
            cache: RwLock::new(HashMap::new()),
        })
    }

    /// Number of structure-bound indices currently retained for reuse.
    #[must_use]
    pub fn cached_index_count(&self) -> usize {
        self.read_cache().len()
    }

    /// Returns this structure snapshot's coordinate slice.
    #[inline]
    fn positions(&self) -> &[[f32; 3]] {
        self.structure.positions()
    }

    /// Resolves one atom position or returns the query out-of-bounds diagnostic.
    #[inline]
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
        let mut selected = Vec::with_capacity(helpers::selection_capacity(universe));

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

            count = count
                .checked_add(1)
                .ok_or_else(|| Diagnostic::new(Code::E3001))?;
        }

        if count == 0 {
            return Err(Diagnostic::new(Code::E5003));
        }

        let divisor = f64::from(count);

        let convert = |value: f64| {
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
    /// A dense nearest-distance buffer keeps pair updates contiguous and avoids
    /// map lookups for atom-indexed nearest-distance state.
    fn iso_layer(
        &self,
        universe: &AtomSelection,
        target: &AtomSelection,
        inner: f32,
        outer: f32,
    ) -> Result<AtomSelection, Diagnostic> {
        validate_shell(inner, outer)?;

        let mut nearest = vec![f32::INFINITY; self.positions().len()];
        let mut reduction_error = None;
        self.for_each_pair(universe, target, outer, |pair| {
            if reduction_error.is_none()
                && let Err(error) = update_nearest_distance(&mut nearest, universe, target, pair)
            {
                reduction_error = Some(error);
            }
        })?;
        if let Some(error) = reduction_error {
            return Err(error);
        }

        select_iso_layer(universe, target, &nearest, inner * inner, outer * outer)
    }

    /// Selects query atoms within `cutoff` of at least one target.
    ///
    /// A dense bit vector records membership directly by atom index and avoids
    /// sorting or deduplicating matched pair endpoints.
    fn within(
        &self,
        query: &AtomSelection,
        target: &AtomSelection,
        cutoff: f32,
    ) -> Result<AtomSelection, Diagnostic> {
        let mut matched = vec![false; self.positions().len()];
        let mut reduction_error = None;
        self.for_each_pair(query, target, cutoff, |pair| {
            if reduction_error.is_none()
                && let Err(error) = mark_pair_match(&mut matched, query, target, pair)
            {
                reduction_error = Some(error);
            }
        })?;
        if let Some(error) = reduction_error {
            return Err(error);
        }

        collect_within(query, target, &matched)
    }

    /// Searches one pair workload while retaining indices in this resolver.
    ///
    /// The resolver itself is configured with `Auto` in a reusable plan.  A
    /// request may still pin a concrete backend; the target index is keyed by
    /// that backend, cutoff and selection, so compatible operations share the
    /// construction without changing the requested algorithm.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic when the requested backend cannot be planned or
    /// the spatial kernel rejects the workload.
    pub fn pairs_with_backend(
        &self,
        query: &AtomSelection,
        target: &AtomSelection,
        cutoff: f32,
        requested_backend: SpatialBackend,
    ) -> Result<Vec<NeighborPair>, Diagnostic> {
        let requested = SpatialSearchOptions {
            backend: requested_backend,
            ..self.options
        };
        let backend =
            cacheable_backend(requested, query.len(), target.len()).map_err(spatial_diagnostic)?;

        match backend {
            SpatialBackend::CellList | SpatialBackend::KdTree => {
                self.cached_pairs(query, target, cutoff, backend)
            }
            _ => {
                let options = SpatialSearchOptions {
                    backend,
                    ..self.options
                };

                pairs_within_with_options(
                    self.positions(),
                    query,
                    target,
                    cutoff,
                    options,
                    self.periodic.as_ref(),
                    self.context,
                )
                .map_err(spatial_diagnostic)
            }
        }
    }

    /// Executes a query using a reusable `CellList` or `KdTree`.
    ///
    /// Atom selections are materialised only where the index APIs require
    /// contiguous integer slices.
    fn cached_pairs(
        &self,
        query: &AtomSelection,
        target: &AtomSelection,
        cutoff: f32,
        backend: SpatialBackend,
    ) -> Result<Vec<NeighborPair>, Diagnostic> {
        let target_indices = collect_selection_indices(target);

        let key = IndexKey {
            backend: backend_key(backend),
            cutoff: index_cutoff_key(backend, cutoff),
            selection: target_indices.into_boxed_slice(),
        };

        let index = self.cached_index(key, backend, cutoff)?;
        let query_indices = collect_selection_indices(query);

        match index.as_ref() {
            CachedIndex::Cell(index) => index
                .pairs(&query_indices, cutoff)
                .map_err(spatial_diagnostic),
            CachedIndex::Kd(index) => index
                .pairs(&query_indices, cutoff)
                .map_err(spatial_diagnostic),
        }
    }

    /// Returns or constructs one structure-bound spatial index.
    ///
    /// Index construction occurs without holding the cache lock. A second
    /// lookup after construction resolves concurrent insertions of the same key
    /// without replacing an index already retained by another thread.
    fn cached_index(
        &self,
        key: IndexKey,
        backend: SpatialBackend,
        cutoff: f32,
    ) -> Result<Arc<CachedIndex<'a>>, Diagnostic> {
        {
            let cache = self.read_cache();

            if let Some(index) = cache.get(&key) {
                return Ok(Arc::clone(index));
            }
        }

        let built = Arc::new(self.build_cached_index(backend, key.selection.as_ref(), cutoff)?);

        let mut cache = self.write_cache();
        let can_retain = cache.len() < MAX_CACHED_INDICES;

        match cache.entry(key) {
            Entry::Occupied(entry) => Ok(Arc::clone(entry.get())),
            Entry::Vacant(entry) if can_retain => {
                entry.insert(Arc::clone(&built));
                Ok(built)
            }
            Entry::Vacant(_) => Ok(built),
        }
    }

    /// Constructs one cacheable index over `targets`.
    ///
    /// K-d trees are independent from the query cutoff, while cell lists retain
    /// the cutoff used to define their spatial grid.
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
                self.periodic,
                self.options.cell_grid,
            )
            .map(CachedIndex::Cell)
            .map_err(spatial_diagnostic),
            SpatialBackend::KdTree => KdTree::build_with_options(
                positions,
                targets,
                self.periodic,
                self.options.kd_periodic,
            )
            .map(CachedIndex::Kd)
            .map_err(spatial_diagnostic),
            _ => Err(Diagnostic::new(Code::E9001)),
        }
    }

    /// Acquires shared access to the spatial-index cache.
    ///
    /// Poisoned locks are recovered because an interrupted cache operation
    /// cannot invalidate the immutable structure snapshot itself.
    fn read_cache(&self) -> RwLockReadGuard<'_, HashMap<IndexKey, Arc<CachedIndex<'a>>>> {
        match self.cache.read() {
            Ok(cache) => cache,
            Err(poisoned) => poisoned.into_inner(),
        }
    }

    /// Acquires exclusive access to mutate the spatial-index cache.
    ///
    /// Poisoned locks are recovered so cache availability is not coupled to a
    /// prior thread panic.
    fn write_cache(&self) -> RwLockWriteGuard<'_, HashMap<IndexKey, Arc<CachedIndex<'a>>>> {
        match self.cache.write() {
            Ok(cache) => cache,
            Err(poisoned) => poisoned.into_inner(),
        }
    }
}

#[cfg(test)]
#[path = "resolver_tests.rs"]
mod tests;
