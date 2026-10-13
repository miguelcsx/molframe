//! Lazy, output-sensitive motif mapping.
//!
//! [`MappingPlan`] holds the per-component candidates and, when the caller asks
//! for *constrained* matching, a cell grid per pairwise `Distance` constraint.
//! [`MappingStream`] then walks the mappings depth first without ever holding
//! more than one tuple, and [`MappingPlan::count`] counts them without
//! materialising any. Enumeration order is identical to the exhaustive
//! [`crate::map_motif`] order restricted to the admitted mappings.
//!
//! [`Admission::Satisfied`] admits a mapping when every `Distance` constraint is
//! intrinsically satisfied, using exactly the arithmetic of the measurement
//! stage, so the admitted set equals the mappings whose measured `satisfied`
//! is `Some(true)` for each distance constraint. [`Admission::Passable`] widens
//! each window to what the verdict profile accepts, which makes the admitted set
//! a superset of every mapping the profile passes.

use crate::mapping::{Candidate, MappingError, component_candidates, referenced_sites};
use crate::numeric::usize_to_u32;
use crate::stream_grid::{Link, intersect};
use crate::{MappedMotif, Motif, VerdictProfile};
use molframe_chem::ComponentProvider;
use molframe_core::contract::AnalysisPolicy;
use molframe_core::execution::CancellationToken;
use molframe_core::index::ResidueIndex;
use molframe_core::structure::Structure;
use std::collections::BTreeMap;
use std::sync::Arc;

/// Nodes visited between cancellation checks.
const CHECK_INTERVAL: u64 = 2048;

/// Mappings of one motif onto one structure, ready to stream or count.
#[derive(Debug)]
pub struct MappingPlan {
    pub(crate) levels: Vec<Box<str>>,
    pub(crate) candidates: Vec<Vec<Candidate>>,
    pub(crate) links: Vec<Vec<Link>>,
}

/// Which mappings a plan admits.
#[derive(Clone, Copy, Debug)]
#[non_exhaustive]
pub enum Admission<'a> {
    /// Every structural mapping, however far apart its atoms lie.
    All,
    /// Mappings whose every distance constraint is within its tolerance.
    Satisfied,
    /// Mappings every distance constraint admits, either within its tolerance
    /// or inside the range the profile's rules on that metric accept, so no
    /// mapping the profile would pass is ever pruned. A constraint no rule
    /// names cannot change a verdict and therefore prunes nothing.
    Passable(&'a VerdictProfile),
}

impl MappingPlan {
    /// Builds candidates and, when `constrained`, the spatial index over the
    /// distance constraints' tolerances ([`Admission::Satisfied`]).
    ///
    /// # Errors
    ///
    /// Returns [`MappingError::Chemistry`] when the provider fails.
    pub fn build(
        structure: &Structure,
        motif: &Motif,
        provider: Option<&dyn ComponentProvider>,
        policy: &AnalysisPolicy,
        constrained: bool,
    ) -> Result<Self, MappingError> {
        let admission = if constrained {
            Admission::Satisfied
        } else {
            Admission::All
        };
        Self::build_admitting(structure, motif, provider, policy, admission)
    }

    /// Builds candidates and the spatial index `admission` calls for.
    ///
    /// # Errors
    ///
    /// Returns [`MappingError::Chemistry`] when the provider fails.
    pub fn build_admitting(
        structure: &Structure,
        motif: &Motif,
        provider: Option<&dyn ComponentProvider>,
        policy: &AnalysisPolicy,
        admission: Admission<'_>,
    ) -> Result<Self, MappingError> {
        let sites = referenced_sites(motif);
        let levels: Vec<Box<str>> = motif.components().keys().cloned().collect();
        let mut candidates = Vec::with_capacity(levels.len());
        for (name, spec) in motif.components() {
            candidates.push(component_candidates(
                structure, name, spec, provider, policy, &sites,
            )?);
        }
        let mut plan = Self {
            links: (0..levels.len()).map(|_| Vec::new()).collect(),
            levels,
            candidates,
        };
        let profile = match admission {
            Admission::All => None,
            Admission::Satisfied => Some(None),
            Admission::Passable(profile) => Some(Some(profile)),
        };
        if let Some(profile) = profile {
            plan.constrain(structure, motif, profile);
        }
        Ok(plan)
    }
}

impl MappingPlan {
    /// Number of motif components, which is the width of every mapping tuple.
    #[must_use]
    pub const fn width(&self) -> usize {
        self.levels.len()
    }

    /// Component names in tuple column order.
    #[must_use]
    pub fn components(&self) -> &[Box<str>] {
        &self.levels
    }

    /// Candidate list of `depth` given the candidates chosen above it.
    ///
    /// `None` means every candidate of the level.
    fn options(&self, depth: usize, chosen: &[u32]) -> Option<Vec<u32>> {
        let mut result: Option<Vec<u32>> = None;
        for link in &self.links[depth] {
            let found = link.reachable(chosen[link.earlier] as usize);
            result = Some(match result {
                Some(previous) => intersect(previous, &found),
                None => found,
            });
        }
        result
    }

    fn residue(&self, depth: usize, candidate: u32) -> ResidueIndex {
        self.candidates[depth][candidate as usize].residue
    }

    fn taken(&self, depth: usize, chosen: &[u32], candidate: u32) -> bool {
        let residue = self.residue(depth, candidate);
        (0..depth).any(|level| self.residue(level, chosen[level]) == residue)
    }

    /// A lazy stream over every admitted mapping, in deterministic order.
    #[must_use]
    pub fn stream(self: &Arc<Self>) -> MappingStream {
        let width = self.levels.len();
        MappingStream {
            plan: Arc::clone(self),
            lists: vec![None; width],
            cursors: vec![0; width],
            chosen: vec![0; width],
            depth: 0,
            started: false,
            done: false,
            steps: 0,
        }
    }

    /// Counts the admitted mappings without materialising any.
    ///
    /// The final component is counted arithmetically, so the cost follows the
    /// number of admitted prefixes, not the number of mappings.
    ///
    /// # Errors
    ///
    /// Returns [`MappingError::Cancelled`] when `cancel` fires.
    pub fn count(&self, cancel: &CancellationToken) -> Result<u64, MappingError> {
        if cancel.is_cancelled() {
            return Err(MappingError::Cancelled);
        }
        if self.levels.is_empty() {
            return Ok(1);
        }
        let mut chosen = vec![0u32; self.levels.len()];
        let mut steps = 0u64;
        self.count_level(0, &mut chosen, &mut steps, cancel)
    }

    fn count_level(
        &self,
        depth: usize,
        chosen: &mut Vec<u32>,
        steps: &mut u64,
        cancel: &CancellationToken,
    ) -> Result<u64, MappingError> {
        *steps += 1;
        if steps.is_multiple_of(CHECK_INTERVAL) && cancel.is_cancelled() {
            return Err(MappingError::Cancelled);
        }
        let options = self.options(depth, chosen);
        let total = self.candidates[depth].len();
        if depth + 1 == self.levels.len() {
            let free = |candidate: u32| !self.taken(depth, chosen, candidate);
            let count = if let Some(list) = &options {
                list.iter().filter(|&&candidate| free(candidate)).count()
            } else {
                {
                    let used = (0..depth)
                        .filter(|&level| {
                            let residue = self.residue(level, chosen[level]);
                            self.candidates[depth]
                                .binary_search_by_key(&residue, |candidate| candidate.residue)
                                .is_ok()
                        })
                        .count();
                    total.saturating_sub(used)
                }
            };
            return Ok(match u64::try_from(count) {
                Ok(count) => count,
                Err(_) => u64::MAX,
            });
        }
        let mut sum = 0u64;
        let length = options.as_ref().map_or(total, Vec::len);
        for position in 0..length {
            let candidate = options
                .as_ref()
                .map_or_else(|| usize_to_u32(position), |list| list[position]);
            if self.taken(depth, chosen, candidate) {
                continue;
            }
            chosen[depth] = candidate;
            sum = sum.saturating_add(self.count_level(depth + 1, chosen, steps, cancel)?);
        }
        Ok(sum)
    }

    /// The residue indices of one mapping tuple, widened to a full mapping.
    #[must_use]
    pub fn mapping(&self, tuple: &[u32]) -> MappedMotif {
        let mut components = BTreeMap::new();
        let mut atoms = BTreeMap::new();
        for (depth, &residue) in tuple.iter().enumerate() {
            let position = match self.candidates[depth]
                .binary_search_by_key(&residue, |candidate| candidate.residue.get())
            {
                Ok(position) => position,
                Err(_) => 0,
            };
            let candidate = &self.candidates[depth][position];
            components.insert(self.levels[depth].clone(), candidate.residue);
            atoms.extend(candidate.atoms.clone());
        }
        MappedMotif { components, atoms }
    }
}

/// A resumable depth-first walk over a [`MappingPlan`].
#[derive(Debug)]
pub struct MappingStream {
    plan: Arc<MappingPlan>,
    lists: Vec<Option<Vec<u32>>>,
    cursors: Vec<usize>,
    chosen: Vec<u32>,
    depth: usize,
    started: bool,
    done: bool,
    steps: u64,
}

impl MappingStream {
    /// Writes the next mapping's residue indices (one per component) into `out`.
    ///
    /// Returns `false` once the stream is exhausted.
    ///
    /// # Errors
    ///
    /// Returns [`MappingError::Cancelled`] when `cancel` fires.
    pub fn next_into(
        &mut self,
        out: &mut [u32],
        cancel: &CancellationToken,
    ) -> Result<bool, MappingError> {
        let width = self.plan.levels.len();
        if self.done {
            return Ok(false);
        }
        if width == 0 {
            self.done = true;
            return Ok(true);
        }
        if self.started {
            self.depth = width - 1;
        } else {
            if cancel.is_cancelled() {
                return Err(MappingError::Cancelled);
            }
            self.started = true;
            self.depth = 0;
            self.lists[0] = None;
            self.cursors[0] = 0;
        }
        loop {
            self.steps += 1;
            if self.steps.is_multiple_of(CHECK_INTERVAL) && cancel.is_cancelled() {
                return Err(MappingError::Cancelled);
            }
            let depth = self.depth;
            let length = self.lists[depth]
                .as_ref()
                .map_or(self.plan.candidates[depth].len(), Vec::len);
            if self.cursors[depth] >= length {
                if depth == 0 {
                    self.done = true;
                    return Ok(false);
                }
                self.depth -= 1;
                continue;
            }
            let position = self.cursors[depth];
            self.cursors[depth] += 1;
            let candidate = match &self.lists[depth] {
                Some(list) => list[position],
                None => match u32::try_from(position) {
                    Ok(position) => position,
                    Err(_) => u32::MAX,
                },
            };
            if self.plan.taken(depth, &self.chosen, candidate) {
                continue;
            }
            self.chosen[depth] = candidate;
            if depth + 1 == width {
                for (level, slot) in out.iter_mut().enumerate() {
                    *slot = self.plan.residue(level, self.chosen[level]).get();
                }
                return Ok(true);
            }
            self.depth += 1;
            self.lists[self.depth] = self.plan.options(self.depth, &self.chosen);
            self.cursors[self.depth] = 0;
        }
    }
}

#[cfg(test)]
#[path = "stream_tests.rs"]
mod tests;
