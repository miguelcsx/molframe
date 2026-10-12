//! `DockQ` and QS-score between structures that do not number atoms alike.
//!
//! The plain scores need the `i`-th atom of the model to answer to the `i`-th
//! of the native. Here chains are matched by sequence identity, residues by
//! alignment inside each chain pair and atoms by name (with chemically
//! equivalent atoms resolved), and the scores run on exactly the atoms that
//! correspond. Chain labels in the arguments are the native structure's.

use crate::CompareError;
use crate::correspondence::{
    ChainAssignment, ChainMapping, assign_chains, map_atoms, map_residues, pair_chains,
};
use crate::docking_quality::{DockQ, DockQOptions, dockq_on_positions};
use crate::interface::{AtomInfo, UNANNOTATED, atom_infos, chain_atoms};
use crate::quaternary::{QsOptions, qs_on_positions};
use crate::workflow::PointMapping;
use molframe_chem::ComponentProvider;
use molframe_core::contract::{Namespace, Status};
use molframe_core::structure::Structure;
use molframe_core::{Code, Diagnostic};
use molframe_seq::Scoring;
use std::collections::BTreeSet;

/// Why a mapped comparison could not be scored.
#[derive(Debug, thiserror::Error)]
pub enum MappedCompareError {
    /// Chains, residues or atoms could not be put in correspondence.
    #[error("correspondence failed: {0}")]
    Mapping(Diagnostic),
    /// The score itself failed on the corresponding atoms.
    #[error(transparent)]
    Compare(#[from] CompareError),
}

/// How the two structures are put in correspondence.
#[derive(Clone, Copy)]
pub struct MappingOptions<'a> {
    /// Chemical component definitions for sequence codes and equivalences.
    pub provider: &'a dyn ComponentProvider,
    /// Namespace in which chain labels are given.
    pub namespace: Namespace,
    /// Alignment scoring for chains and residues.
    pub scoring: Scoring,
    /// Smallest sequence identity at which two chains correspond, in `[0, 1]`.
    pub min_identity: f64,
    /// Bound on the chemically equivalent atom mappings tried per residue.
    pub automorphism_limit: usize,
}

impl std::fmt::Debug for MappingOptions<'_> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("MappingOptions")
            .field("provider", &self.provider.version())
            .field("namespace", &self.namespace)
            .field("scoring", &self.scoring)
            .field("min_identity", &self.min_identity)
            .field("automorphism_limit", &self.automorphism_limit)
            .finish()
    }
}

/// What corresponded for the scored chains.
#[derive(Clone, Debug, PartialEq)]
pub struct MappedComparison {
    /// Native chain to model chain, for the chains that were scored.
    pub chains: Vec<ChainMapping>,
    /// Native atom to model atom, in native residue order, for those chains.
    pub atoms: PointMapping,
    /// Residue pairs whose equivalent atoms were swapped to fit better.
    pub swapped_residues: usize,
    /// Result status; tied atom correspondences are ambiguous.
    pub status: Status,
    /// Minimum-distance alternatives for chemically equivalent atom mappings.
    pub atom_alternatives: Vec<crate::correspondence::AtomMappingAlternative>,
    /// How many chain assignments of equal sequence identity were scored. More
    /// than one means the sequences alone could not say which model chain
    /// answers to which native chain (a homo-oligomer), and the reported
    /// assignment is the one that scored best.
    pub chain_assignments_tried: usize,
}

/// The corresponding rows of both structures, and which rows belong to each
/// requested native chain.
struct Rows {
    model: Vec<[f32; 3]>,
    native: Vec<[f32; 3]>,
    chains: Vec<Vec<usize>>,
    /// Annotation of each row, taken from the native atom it came from.
    info: Vec<AtomInfo>,
}

/// Targets a native chain may answer to: the one the optimal assignment chose
/// first, then every other model chain of exactly its sequence identity.
fn tied_targets(assignment: &ChainAssignment, label: &str) -> Vec<ChainMapping> {
    let Some(primary) = assignment
        .primary
        .iter()
        .find(|mapping| mapping.reference == label)
    else {
        return unassigned_targets(assignment, label);
    };
    let mut targets = vec![primary.clone()];
    targets.extend(
        assignment
            .alternatives
            .iter()
            .filter(|alternative| {
                alternative.reference == label
                    && alternative.identity.total_cmp(&primary.identity).is_eq()
            })
            .map(|alternative| ChainMapping {
                reference: alternative.reference.clone(),
                target: alternative.target.clone(),
                identity: alternative.identity,
            }),
    );
    targets
}

/// Targets for a native chain the one-to-one assignment left out, which happens
/// when the native holds more identical chains than the model. Every model
/// chain of the best identity this chain reached is a candidate; the scores
/// then say which one the native chain really answers to.
fn unassigned_targets(assignment: &ChainAssignment, label: &str) -> Vec<ChainMapping> {
    let own = || {
        assignment
            .alternatives
            .iter()
            .filter(|alternative| alternative.reference == label)
    };
    let Some(best) = own()
        .map(|alternative| alternative.identity)
        .max_by(f64::total_cmp)
    else {
        return Vec::new();
    };
    own()
        .filter(|alternative| alternative.identity.total_cmp(&best).is_eq())
        .map(|alternative| ChainMapping {
            reference: alternative.reference.clone(),
            target: alternative.target.clone(),
            identity: alternative.identity,
        })
        .collect()
}

/// Every assignment of the two native chains to distinct tied model chains,
/// the sequence-optimal one first.
fn candidate_assignments(
    assignment: &ChainAssignment,
    labels: [&str; 2],
) -> Vec<[ChainMapping; 2]> {
    let first = tied_targets(assignment, labels[0]);
    let second = tied_targets(assignment, labels[1]);
    let mut candidates = Vec::with_capacity(first.len() * second.len());
    for left in &first {
        for right in &second {
            if labels[0] == labels[1] || left.target != right.target {
                candidates.push([left.clone(), right.clone()]);
            }
        }
    }
    candidates
}

fn correspond(
    model: &Structure,
    native: &Structure,
    mapping: &MappingOptions<'_>,
    chains: Vec<ChainMapping>,
    chain_assignments_tried: usize,
) -> Result<MappedComparison, MappedCompareError> {
    let residues = map_residues(
        native,
        model,
        &chains,
        mapping.provider,
        mapping.namespace,
        mapping.scoring,
    )
    .map_err(MappedCompareError::Mapping)?;
    let atoms = map_atoms(
        native,
        model,
        &residues,
        mapping.provider,
        mapping.automorphism_limit,
    )?;
    Ok(MappedComparison {
        chains,
        atoms: atoms.mapping,
        swapped_residues: atoms.swapped_residues,
        status: atoms.status,
        atom_alternatives: atoms.alternatives,
        chain_assignments_tried,
    })
}

fn rows(
    model: &Structure,
    native: &Structure,
    comparison: &MappedComparison,
    labels: &[&str],
    namespace: Namespace,
) -> Result<Rows, MappedCompareError> {
    let wanted: Vec<BTreeSet<usize>> = labels
        .iter()
        .map(|label| chain_atoms(native, label, namespace).map(|atoms| atoms.into_iter().collect()))
        .collect::<Result<_, _>>()?;
    let mut rows = Rows {
        model: Vec::new(),
        native: Vec::new(),
        chains: vec![Vec::new(); labels.len()],
        info: Vec::new(),
    };
    let infos = atom_infos(native);
    for pair in comparison.atoms.matches() {
        let (Some(from), Some(to)) = (
            native.positions().get(pair.reference),
            model.positions().get(pair.model),
        ) else {
            continue;
        };
        let row = rows.native.len();
        rows.native.push(*from);
        rows.model.push(*to);
        rows.info.push(match infos.get(pair.reference) {
            Some(&info) => info,
            None => UNANNOTATED,
        });
        for (chain, atoms) in rows.chains.iter_mut().zip(&wanted) {
            if atoms.contains(&pair.reference) {
                chain.push(row);
            }
        }
    }
    if rows.chains.iter().any(Vec::is_empty) {
        return Err(MappedCompareError::Compare(CompareError::NoComparablePairs));
    }
    Ok(rows)
}

/// Scores every candidate chain assignment of `labels` and keeps the best.
///
/// Both scores are independent of the frame the structures share, so equal
/// sequences are told apart by how well the resulting correspondence scores,
/// never by where the chains sit. Cost: one residue and atom correspondence and
/// one score per candidate; `k` identical chains give `k * (k - 1)` candidates.
/// The sequence-optimal assignment is tried first and wins ties.
fn best_assignment<S>(
    model: &Structure,
    native: &Structure,
    mapping: &MappingOptions<'_>,
    labels: [&str; 2],
    score: impl Fn(&Rows) -> Result<S, CompareError>,
    value: impl Fn(&S) -> f64,
) -> Result<(S, MappedComparison), MappedCompareError> {
    let assignment = assign_chains(
        native,
        model,
        mapping.provider,
        mapping.namespace,
        mapping.scoring,
        mapping.min_identity,
    )
    .map_err(MappedCompareError::Mapping)?;
    let candidates = candidate_assignments(&assignment, labels);
    let tried = candidates.len();
    let mut best: Option<(S, MappedComparison)> = None;
    let mut first_error = None;
    for chains in candidates {
        let attempt =
            correspond(model, native, mapping, chains.into(), tried).and_then(|comparison| {
                let rows = rows(model, native, &comparison, &labels, mapping.namespace)?;
                Ok((score(&rows)?, comparison))
            });
        match attempt {
            Ok((scored, comparison)) => {
                let better = best
                    .as_ref()
                    .is_none_or(|(held, _)| value(&scored) > value(held));
                if better {
                    best = Some((scored, comparison));
                }
            }
            Err(error) => {
                first_error.get_or_insert(error);
            }
        }
    }
    match (best, first_error) {
        (Some(found), _) => Ok(found),
        (None, Some(error)) => Err(error),
        (None, None) => Err(MappedCompareError::Compare(CompareError::NoComparablePairs)),
    }
}

/// Scores a docking model against its native complex across renamed and
/// renumbered chains.
///
/// `receptor` and `ligand` are native chain labels. Where several model chains
/// match a native chain equally well by sequence, the assignment that scores
/// best is used and [`MappedComparison::chain_assignments_tried`] says how many
/// were compared.
///
/// # Errors
///
/// Returns a mapping failure when chains, residues or atoms cannot be matched,
/// and the errors of [`dockq`](crate::dockq) for the score itself. A named
/// chain with no corresponding atoms is [`CompareError::NoComparablePairs`].
pub fn mapped_dockq(
    model: &Structure,
    native: &Structure,
    receptor: &str,
    ligand: &str,
    mapping: &MappingOptions<'_>,
    options: DockQOptions,
) -> Result<(DockQ, MappedComparison), MappedCompareError> {
    best_assignment(
        model,
        native,
        mapping,
        [receptor, ligand],
        |rows| {
            dockq_on_positions(
                &rows.model,
                &rows.native,
                &rows.chains[0],
                &rows.chains[1],
                &rows.info,
                options,
            )
        },
        |score| score.score,
    )
}

/// Scores a docking model with the model chains named by the caller.
///
/// `receptor` and `ligand` are `(native chain, model chain)` pairs. No chain
/// assignment is searched: the pairing is taken as given, so a pairing the
/// sequences would not have chosen is scored rather than refused, and each
/// pair's [`ChainMapping::identity`] reports how well its sequences agree.
/// Pairs below [`MappingOptions::min_identity`] are refused. Because nothing is
/// assigned one to one, this also works when the native holds more identical
/// chains than the model.
///
/// # Errors
///
/// As [`mapped_dockq`], and a mapping failure when a named chain is missing or a
/// pair is below the identity threshold.
pub fn mapped_dockq_for_chains(
    model: &Structure,
    native: &Structure,
    receptor: (&str, &str),
    ligand: (&str, &str),
    mapping: &MappingOptions<'_>,
    options: DockQOptions,
) -> Result<(DockQ, MappedComparison), MappedCompareError> {
    let pair = |labels: (&str, &str)| {
        let found = pair_chains(
            native,
            model,
            labels,
            mapping.provider,
            mapping.namespace,
            mapping.scoring,
        )
        .map_err(MappedCompareError::Mapping)?;
        if found.identity < mapping.min_identity {
            return Err(MappedCompareError::Mapping(
                Diagnostic::new(Code::E4003)
                    .with_context("chain", labels.0)
                    .with_context("identity", found.identity.to_string())
                    .with_context("min_identity", mapping.min_identity.to_string()),
            ));
        }
        Ok(found)
    };
    let chains = vec![pair(receptor)?, pair(ligand)?];
    let comparison = correspond(model, native, mapping, chains, 1)?;
    let rows = rows(
        model,
        native,
        &comparison,
        &[receptor.0, ligand.0],
        mapping.namespace,
    )?;
    let scored = dockq_on_positions(
        &rows.model,
        &rows.native,
        &rows.chains[0],
        &rows.chains[1],
        &rows.info,
        options,
    )?;
    Ok((scored, comparison))
}

/// Scores the contact overlap of one chain pair across renamed and
/// renumbered chains.
///
/// `first` and `second` are native chain labels; tied chains are resolved as in
/// [`mapped_dockq`].
///
/// # Errors
///
/// As [`mapped_dockq`], with the errors of [`qs_score`](crate::qs_score) for
/// the score.
pub fn mapped_qs_score(
    model: &Structure,
    native: &Structure,
    first: &str,
    second: &str,
    mapping: &MappingOptions<'_>,
    options: QsOptions,
) -> Result<(f64, MappedComparison), MappedCompareError> {
    best_assignment(
        model,
        native,
        mapping,
        [first, second],
        |rows| {
            qs_on_positions(
                &rows.model,
                &rows.native,
                &rows.chains[0],
                &rows.chains[1],
                &rows.info,
                options,
            )
        },
        |score| *score,
    )
}

#[cfg(test)]
#[path = "mapped_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "mapped_pairing_tests.rs"]
mod pairing_tests;
