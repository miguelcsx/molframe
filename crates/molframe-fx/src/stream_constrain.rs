//! Spatial pruning of a mapping plan by its pairwise distance constraints.
//!
//! Each pairwise `Distance` constraint becomes a cell grid over the later
//! component's candidate positions, so enumeration visits only the partners
//! within reach. Constraints between two sites of one component filter that
//! component's candidates directly.

use crate::Comparison;
use crate::stream::MappingPlan;
use crate::stream_grid::{Grid, Link, Window, site_points};
use crate::{AtomSite, Constraint, Motif, NamedConstraint, VerdictProfile, VerdictRule};
use molframe_core::structure::Structure;

/// The closed range of distances a rule can accept.
fn accepted_by(comparison: Comparison) -> (f64, f64) {
    match comparison {
        Comparison::LessThan(limit) | Comparison::AtMost(limit) => (f64::NEG_INFINITY, limit),
        Comparison::GreaterThan(limit) | Comparison::AtLeast(limit) => (limit, f64::INFINITY),
        Comparison::Between(low, high) => (low, high),
    }
}

/// Distances a profile can pass for `metric`: every rule on it must accept.
///
/// A metric no rule mentions cannot change a verdict, so every distance is
/// passable and the constraint prunes nothing.
fn passable(profile: &VerdictProfile, metric: &str) -> (f64, f64) {
    profile
        .rules()
        .iter()
        .filter(|rule: &&VerdictRule| rule.metric.as_ref() == metric)
        .map(|rule| accepted_by(rule.comparison))
        .fold(
            (f64::NEG_INFINITY, f64::INFINITY),
            |(low, high), (other_low, other_high)| (low.max(other_low), high.min(other_high)),
        )
}

fn window_of(named: &NamedConstraint, profile: Option<&VerdictProfile>) -> Option<Window> {
    let Constraint::Distance {
        target, tolerance, ..
    } = &named.constraint
    else {
        return None;
    };
    Some(Window {
        target: *target,
        tolerance: *tolerance,
        accepted: profile.map(|profile| passable(profile, &named.name)),
    })
}

impl MappingPlan {
    fn level_of(&self, name: &str) -> Option<usize> {
        self.levels.iter().position(|level| level.as_ref() == name)
    }

    pub(crate) fn constrain(
        &mut self,
        structure: &Structure,
        motif: &Motif,
        profile: Option<&VerdictProfile>,
    ) {
        let pairs: Vec<_> = motif
            .constraints()
            .iter()
            .filter_map(|named| {
                let window = window_of(named, profile)?;
                if window.admits_everything() {
                    return None;
                }
                let Constraint::Distance { first, second, .. } = &named.constraint else {
                    return None;
                };
                let a = self.level_of(&first.component)?;
                let b = self.level_of(&second.component)?;
                Some((a, b, first, second, window))
            })
            .collect();
        for &(a, b, first, second, window) in &pairs {
            if a == b {
                self.filter_internal(structure, a, first, second, window);
            }
        }
        for &(a, b, first, second, window) in &pairs {
            if a == b {
                continue;
            }
            let (earlier, earlier_site, later, later_site) = if a < b {
                (a, first, b, second)
            } else {
                (b, second, a, first)
            };
            self.add_link(
                structure,
                (earlier, earlier_site),
                (later, later_site),
                window,
            );
        }
    }

    fn filter_internal(
        &mut self,
        structure: &Structure,
        level: usize,
        first: &AtomSite,
        second: &AtomSite,
        window: Window,
    ) {
        self.candidates[level].retain(|candidate| {
            let left = site_points(structure, candidate, first);
            let right = site_points(structure, candidate, second);
            left.iter().any(|(atom_a, p)| {
                right.iter().any(|(atom_b, q)| {
                    (atom_a != atom_b || first == second)
                        && window.admits(molframe_geom::distance(*p, *q))
                })
            })
        });
    }

    fn add_link(
        &mut self,
        structure: &Structure,
        (earlier, earlier_site): (usize, &AtomSite),
        (later, later_site): (usize, &AtomSite),
        window: Window,
    ) {
        let reach = window.reach();
        let cell = if reach.is_finite() {
            (reach.max(0.0) * (1.0 + 1e-6) + 1e-6).max(1e-6)
        } else {
            f64::INFINITY
        };
        let grid = Grid::build(
            cell,
            self.candidates[later]
                .iter()
                .enumerate()
                .flat_map(|(index, candidate)| {
                    let index = u32::try_from(index).unwrap_or(u32::MAX);
                    site_points(structure, candidate, later_site)
                        .into_iter()
                        .map(move |(_, point)| (index, point))
                }),
        );
        let earlier_points = self.candidates[earlier]
            .iter()
            .map(|candidate| site_points(structure, candidate, earlier_site))
            .collect();
        self.links[later].push(Link {
            earlier,
            earlier_points,
            window,
            grid,
        });
    }
}
