//! Cell grid over candidate site positions for pairwise distance pruning.

use crate::AtomSite;
use crate::mapping::Candidate;
use molframe_core::index::AtomIndex;
use molframe_core::structure::Structure;
use std::collections::HashMap;

pub(crate) type Point = [f32; 3];
pub(crate) type Cell = [i64; 3];

pub(crate) fn satisfied(distance: f64, target: f64, tolerance: f64) -> bool {
    (distance - target).abs() <= tolerance
}

/// Positions of the alternatives of one site on one candidate.
pub(crate) type SitePoints = Vec<(AtomIndex, Point)>;

#[derive(Debug)]
pub(crate) struct Grid {
    pub(crate) cell: f64,
    pub(crate) cells: HashMap<Cell, Vec<(u32, Point)>>,
}

impl Grid {
    pub(crate) fn key(&self, point: Point) -> Cell {
        let scaled = |value: f32| {
            let cell = (f64::from(value) / self.cell).floor();
            // Saturating cast: a coordinate this far out cannot be near anything.
            #[allow(clippy::cast_possible_truncation)]
            let cell = cell as i64;
            cell
        };
        [scaled(point[0]), scaled(point[1]), scaled(point[2])]
    }
}

/// One pairwise distance constraint seen from the later component.
#[derive(Debug)]
pub(crate) struct Link {
    pub(crate) earlier: usize,
    /// Points of the earlier site, per earlier candidate.
    pub(crate) earlier_points: Vec<SitePoints>,
    pub(crate) target: f64,
    pub(crate) tolerance: f64,
    pub(crate) grid: Grid,
}

pub(crate) fn site_points(
    structure: &Structure,
    candidate: &Candidate,
    site: &AtomSite,
) -> SitePoints {
    match candidate.atoms.get(site) {
        Some(atoms) => atoms
            .iter()
            .filter_map(|&atom| {
                structure
                    .data()
                    .atom(atom)
                    .and_then(molframe_core::structure::AtomRef::position)
                    .map(|position| (atom, molframe_geom::Rigid::IDENTITY.apply(position)))
            })
            .collect(),
        None => SitePoints::new(),
    }
}

impl Link {
    /// Later-level candidate indices within the constraint of `earlier_candidate`, sorted.
    pub(crate) fn reachable(&self, earlier_candidate: usize) -> Vec<u32> {
        let mut found = Vec::new();
        for (_, point) in &self.earlier_points[earlier_candidate] {
            let centre = self.grid.key(*point);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let key = [centre[0] + dx, centre[1] + dy, centre[2] + dz];
                        let Some(entries) = self.grid.cells.get(&key) else {
                            continue;
                        };
                        for (index, other) in entries {
                            if satisfied(
                                molframe_geom::distance(*point, *other),
                                self.target,
                                self.tolerance,
                            ) {
                                found.push(*index);
                            }
                        }
                    }
                }
            }
        }
        found.sort_unstable();
        found.dedup();
        found
    }
}

pub(crate) fn intersect(left: Vec<u32>, right: &[u32]) -> Vec<u32> {
    left.into_iter()
        .filter(|value| right.binary_search(value).is_ok())
        .collect()
}
