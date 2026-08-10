//! Membrane leaflet identification from an explicit lipid-site graph.

use pdbiox_core::selection::AtomSelection;
use pdbiox_spatial::{PeriodicBox, SpatialBackend, SpatialError, pairs_within};
use std::collections::BTreeMap;

/// Explicit graph construction policy for leaflet identification.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LeafletOptions {
    /// Maximum separation between connected lipid representative sites.
    pub connection_distance: f32,
    /// Spatial implementation to use.
    pub backend: SpatialBackend,
}

/// One connected membrane leaflet in deterministic atom-index order.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Leaflet {
    /// Representative site indices assigned to this component.
    pub sites: Vec<u32>,
}

/// Identifies connected components of caller-selected lipid representative sites.
///
/// Component membership is purely geometric. Lipid recognition and the choice
/// of a representative head-group site remain explicit chemistry/query steps,
/// so this kernel contains no residue-name or atom-name tables.
///
/// # Errors
///
/// Returns [`SpatialError`] for an invalid connection distance, out-of-range
/// representative index, or invalid periodic/spatial input.
pub fn identify_leaflets(
    positions: &[[f32; 3]],
    lipid_sites: &AtomSelection,
    options: LeafletOptions,
    periodic: Option<&PeriodicBox>,
) -> Result<Vec<Leaflet>, SpatialError> {
    let pairs = pairs_within(
        positions,
        lipid_sites,
        lipid_sites,
        options.connection_distance,
        options.backend,
        periodic,
    )?;
    let sites: Vec<u32> = lipid_sites.into_iter().collect();
    let mut disjoint = DisjointSet::new(sites.len());
    let lookup: BTreeMap<u32, usize> = sites
        .iter()
        .copied()
        .enumerate()
        .map(|(local, atom)| (atom, local))
        .collect();

    for pair in pairs {
        let Some(&left) = lookup.get(&pair.first) else {
            continue;
        };
        let Some(&right) = lookup.get(&pair.second) else {
            continue;
        };
        disjoint.join(left, right);
    }

    let mut components: BTreeMap<usize, Vec<u32>> = BTreeMap::new();
    for (local, atom) in sites.into_iter().enumerate() {
        components
            .entry(disjoint.root(local))
            .or_default()
            .push(atom);
    }
    let mut leaflets: Vec<Leaflet> = components
        .into_values()
        .map(|sites| Leaflet { sites })
        .collect();
    leaflets.sort_by_key(|leaflet| leaflet.sites.first().copied());
    Ok(leaflets)
}

struct DisjointSet {
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl DisjointSet {
    fn new(size: usize) -> Self {
        Self {
            parent: (0..size).collect(),
            rank: vec![0; size],
        }
    }

    fn root(&mut self, item: usize) -> usize {
        let parent = self.parent[item];
        if parent != item {
            self.parent[item] = self.root(parent);
        }
        self.parent[item]
    }

    fn join(&mut self, left: usize, right: usize) {
        let mut left_root = self.root(left);
        let mut right_root = self.root(right);
        if left_root == right_root {
            return;
        }
        if self.rank[left_root] < self.rank[right_root] {
            std::mem::swap(&mut left_root, &mut right_root);
        }
        self.parent[right_root] = left_root;
        if self.rank[left_root] == self.rank[right_root] {
            self.rank[left_root] += 1;
        }
    }
}

#[cfg(test)]
#[path = "leaflet_tests.rs"]
mod tests;
