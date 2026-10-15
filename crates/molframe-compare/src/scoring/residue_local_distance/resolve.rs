//! Choosing the naming of interchangeable atoms that preserves most distances.
//!
//! Every residue with an admitted exchange is a site. A site's choice is judged
//! only by the reference distances from the atoms the exchange moves to atoms of
//! other residues, which is all the exchange can change. Sites are first judged
//! at once against the model's own naming, which cannot depend on visiting
//! order, and then revisited one at a time against each other's current choice
//! until a sweep changes nothing or the sweep bound is reached. A change needs a
//! strictly better count, so ties keep the earlier choice and the outcome is
//! deterministic. Neighbour lists are kept only for moved atoms.

use super::options::ResidueLddtOptions;
use super::score::{eligible_pairs, preserved};
use super::symmetry::{Permutations, SymmetryCache};
use super::table::Table;
use crate::mapped::MappedCompareError;
use molframe_chem::ComponentProvider;
use molframe_core::ExecutionContext;
use std::collections::BTreeMap;

const SWEEPS: usize = 8;

pub(super) struct Outcome {
    /// For each reference atom, the atom whose model position it is given.
    pub assigned: Vec<u32>,
    pub renamed: usize,
    pub unresolved: usize,
}

struct Site {
    fragment: Vec<u32>,
    moved: Vec<(u32, usize)>,
    permutations: Permutations,
    choice: usize,
}

pub(super) fn resolve(
    table: &Table<'_>,
    provider: &dyn ComponentProvider,
    options: &ResidueLddtOptions,
    limit: usize,
    context: &ExecutionContext,
) -> Result<Outcome, MappedCompareError> {
    let mut assigned: Vec<u32> = table.identity();
    let (mut sites, unresolved) = sites_of(table, provider, limit)?;
    if sites.is_empty() {
        return Ok(Outcome {
            assigned,
            renamed: 0,
            unresolved,
        });
    }
    let neighbours = neighbours_of(table, &sites, options, context)?;

    let baseline = assigned.clone();
    let choices: Vec<usize> = sites
        .iter()
        .map(|site| best(site, &neighbours, table, &baseline, options, 0))
        .collect();
    for (site, choice) in sites.iter_mut().zip(choices) {
        apply(site, choice, &mut assigned);
    }
    for _ in 0..SWEEPS {
        let mut changed = false;
        for site in &mut sites {
            let current = site.choice;
            let choice = best(site, &neighbours, table, &assigned, options, current);
            if choice != current {
                apply(site, choice, &mut assigned);
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    let renamed = sites.iter().filter(|site| site.choice != 0).count();
    Ok(Outcome {
        assigned,
        renamed,
        unresolved,
    })
}

fn sites_of(
    table: &Table<'_>,
    provider: &dyn ComponentProvider,
    limit: usize,
) -> Result<(Vec<Site>, usize), MappedCompareError> {
    let mut cache = SymmetryCache::default();
    let mut components = BTreeMap::new();
    let mut sites = Vec::new();
    let mut unresolved = 0;
    for (entry, span) in table.spans.iter().enumerate() {
        let name = table.components[entry];
        if name.is_empty() {
            continue;
        }
        if !components.contains_key(name) {
            let found = provider.get(name).map_err(MappedCompareError::Mapping)?;
            components.insert(name, found);
        }
        let Some(Some(component)) = components.get(name) else {
            continue;
        };
        let mut present = Vec::new();
        let mut fragment = Vec::new();
        for (position, atom) in component.atoms.iter().enumerate() {
            let found =
                (span.start..span.end).find(|row| table.names[*row as usize] == &*atom.name);
            if let (Some(row), Ok(position)) = (found, u32::try_from(position)) {
                present.push(position);
                fragment.push(row);
            }
        }
        if fragment.len() < 2 {
            continue;
        }
        let Some(permutations) = cache.permutations(component, &present, limit) else {
            unresolved += 1;
            continue;
        };
        let mut moved: Vec<(u32, usize)> = Vec::new();
        for permutation in permutations.iter() {
            for (position, target) in permutation.iter().enumerate() {
                if *target as usize != position && !moved.iter().any(|(_, seen)| *seen == position)
                {
                    moved.push((fragment[position], position));
                }
            }
        }
        if !moved.is_empty() {
            moved.sort_unstable();
            sites.push(Site {
                fragment,
                moved,
                permutations,
                choice: 0,
            });
        }
    }
    Ok((sites, unresolved))
}

/// Reference distances from moved atoms to atoms of other residues.
fn neighbours_of(
    table: &Table<'_>,
    sites: &[Site],
    options: &ResidueLddtOptions,
    context: &ExecutionContext,
) -> Result<BTreeMap<u32, Vec<(u32, f32)>>, MappedCompareError> {
    let mut moved = vec![false; table.reference.len()];
    for site in sites {
        for (atom, _) in &site.moved {
            moved[*atom as usize] = true;
        }
    }
    let rows = eligible_pairs(table, options, context, |first, second| {
        moved[first as usize] || moved[second as usize]
    })?;
    let mut neighbours: BTreeMap<u32, Vec<(u32, f32)>> = BTreeMap::new();
    for (first, second, distance) in rows {
        if moved[first as usize] {
            neighbours
                .entry(first)
                .or_default()
                .push((second, distance));
        }
        if moved[second as usize] {
            neighbours
                .entry(second)
                .or_default()
                .push((first, distance));
        }
    }
    Ok(neighbours)
}

fn target(site: &Site, choice: usize, atom: u32, position: usize) -> u32 {
    match choice {
        0 => atom,
        choice => site.fragment[site.permutations[choice - 1][position] as usize],
    }
}

/// The choice with the most preserved tests, keeping `current` unless beaten.
fn best(
    site: &Site,
    neighbours: &BTreeMap<u32, Vec<(u32, f32)>>,
    table: &Table<'_>,
    assigned: &[u32],
    options: &ResidueLddtOptions,
    current: usize,
) -> usize {
    let score = |choice: usize| -> u64 {
        let mut total = 0;
        for (atom, position) in &site.moved {
            let Some(own) = table.model[target(site, choice, *atom, *position) as usize] else {
                continue;
            };
            let Some(rows) = neighbours.get(atom) else {
                continue;
            };
            for (other, distance) in rows {
                if let Some(partner) = table.model[assigned[*other as usize] as usize] {
                    total += preserved(&options.tolerances, f64::from(*distance), own, partner);
                }
            }
        }
        total
    };
    let mut chosen = current;
    let mut top = score(current);
    for choice in 0..=site.permutations.len() {
        let candidate = score(choice);
        if candidate > top {
            top = candidate;
            chosen = choice;
        }
    }
    chosen
}

fn apply(site: &mut Site, choice: usize, assigned: &mut [u32]) {
    for (position, atom) in site.fragment.iter().enumerate() {
        assigned[*atom as usize] = target(site, choice, *atom, position);
    }
    site.choice = choice;
}
