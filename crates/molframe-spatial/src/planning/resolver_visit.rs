//! Streaming pair dispatch for structure-bound spatial reductions.

use super::helpers::{
    backend_key, cacheable_backend, collect_selection_indices, index_cutoff_key, spatial_diagnostic,
};
use super::{CachedIndex, IndexKey, StructureSpatial};
use crate::{
    NeighborPair, PairQuery, SpatialBackend, SpatialSearchOptions, for_each_pairs_within_unsorted,
};
use molframe_core::diagnostic::Diagnostic;
use molframe_core::selection::AtomSelection;

impl StructureSpatial<'_> {
    /// Streams one pair workload into a reducer when the backend supports it.
    ///
    /// Cell-list and brute-force searches never retain matches. Backends whose
    /// public index still returns a vector are forwarded without changing their
    /// established planning or error behaviour.
    pub(super) fn for_each_pair(
        &self,
        query: &AtomSelection,
        target: &AtomSelection,
        cutoff: f32,
        mut emit: impl FnMut(NeighborPair),
    ) -> Result<(), Diagnostic> {
        if self.periodic.is_some() {
            return for_each_pairs_within_unsorted(
                &PairQuery {
                    positions: self.positions(),
                    left: query,
                    right: target,
                    cutoff,
                    options: self.options,
                    periodic: self.periodic.as_ref(),
                    context: self.context,
                },
                emit,
            )
            .map_err(spatial_diagnostic);
        }

        let backend = cacheable_backend(self.options, query.len(), target.len())
            .map_err(spatial_diagnostic)?;
        match backend {
            SpatialBackend::CellList | SpatialBackend::KdTree => {
                self.cached_for_each_pair(query, target, cutoff, backend, &mut emit)
            }
            _ => for_each_pairs_within_unsorted(
                &PairQuery {
                    positions: self.positions(),
                    left: query,
                    right: target,
                    cutoff,
                    options: SpatialSearchOptions {
                        backend,
                        ..self.options
                    },
                    periodic: None,
                    context: self.context,
                },
                emit,
            )
            .map_err(spatial_diagnostic),
        }
    }

    /// Streams cached-index matches into a query reduction.
    fn cached_for_each_pair(
        &self,
        query: &AtomSelection,
        target: &AtomSelection,
        cutoff: f32,
        backend: SpatialBackend,
        mut emit: impl FnMut(NeighborPair),
    ) -> Result<(), Diagnostic> {
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
                .for_each_pair(&query_indices, cutoff, emit)
                .map_err(spatial_diagnostic),
            CachedIndex::Kd(index) => {
                for pair in index
                    .pairs(&query_indices, cutoff)
                    .map_err(spatial_diagnostic)?
                {
                    emit(pair);
                }
                Ok(())
            }
        }
    }
}
