//! Residues that form the interface between two chains.
//!
//! An interface residue is one with any atom within the cutoff of an atom in the
//! other chain. The residues of both chains that meet this test are returned
//! together, sorted, so the result names the whole contact patch rather than one
//! side of it.
//!
//! The atom pairs come from the shared spatial search restricted to the two
//! chains, so the cost tracks the size of the smaller chain's neighbourhood.

use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::selection::AtomSelection;
use molframe_core::structure::{AtomRef, ChainRef, Structure};
use molframe_core::{ExecutionContext, index::ResidueIndex};
use molframe_spatial::{
    PairQuery, PeriodicBox, SpatialBackend, SpatialError, SpatialSearchOptions, StructureSpatial,
    reduce_pairs_within_unsorted,
};

/// Which chain identifier namespace a chain name is read in.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ChainNamespace {
    /// The normalised label identifier (`label_asym_id`).
    Label,
    /// The depositor's identifier (`auth_asym_id`).
    Auth,
}

/// A chain name and, optionally, the namespace it is written in.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ChainSelector<'a> {
    /// The identifier text.
    pub name: &'a str,
    /// The namespace to match in; `None` accepts either, provided the name
    /// means the same chains in both.
    pub namespace: Option<ChainNamespace>,
}

impl<'a> ChainSelector<'a> {
    /// A name matched in either namespace, rejected when the two disagree.
    #[must_use]
    pub const fn any(name: &'a str) -> Self {
        Self {
            name,
            namespace: None,
        }
    }

    /// A name matched only against label identifiers.
    #[must_use]
    pub const fn label(name: &'a str) -> Self {
        Self {
            name,
            namespace: Some(ChainNamespace::Label),
        }
    }

    /// A name matched only against author identifiers.
    #[must_use]
    pub const fn auth(name: &'a str) -> Self {
        Self {
            name,
            namespace: Some(ChainNamespace::Auth),
        }
    }
}

/// Inputs of an interface search.
#[derive(Clone, Copy, Debug)]
pub struct InterfaceOptions<'a> {
    /// The first chain.
    pub first: ChainSelector<'a>,
    /// The second chain.
    pub second: ChainSelector<'a>,
    /// Contact cutoff in ångström.
    pub cutoff: f32,
    /// Spatial implementation to use.
    pub backend: SpatialBackend,
    /// Apply the structure unit cell and minimum-image convention.
    pub periodic: bool,
}

/// Why an interface could not be evaluated.
#[derive(Debug, thiserror::Error)]
pub enum InterfaceError {
    /// The neighbour search failed.
    #[error(transparent)]
    Spatial(#[from] SpatialError),
    /// A name written without a namespace denotes different chains as a label
    /// and as an author identifier, so it cannot be resolved.
    #[error(
        "chain name {0:?} is ambiguous: it names different chains as a label and as an author id; choose a namespace"
    )]
    AmbiguousChain(String),
    /// Both selectors resolve to the same chain, which has no interface with itself.
    #[error("the two chain selectors resolve to the same chain")]
    SameChain,
    /// Periodic geometry was requested without a unit cell.
    #[error("periodic interface search requires a unit cell")]
    MissingCell,
    /// The structure carries only the placeholder cell some files write when
    /// they have none, which would invent periodic images.
    #[error("the unit cell is a placeholder, not a real crystallographic cell")]
    PlaceholderCell,
}

impl From<InterfaceError> for Diagnostic {
    fn from(error: InterfaceError) -> Self {
        match error {
            InterfaceError::Spatial(error) => error.into(),
            other => Self::new(Code::E4002)
                .with_context("interface", "invalid_selection")
                .with_message(other.to_string()),
        }
    }
}

/// Returns the residues at the interface between the two named chains.
///
/// A name is matched in either the label or the author namespace. When it
/// denotes different chains in the two namespaces the call fails rather than
/// guessing; use [`chain_interface_with_options`] to pick a namespace. When
/// either chain is absent or the two never approach within the cutoff, the
/// result is empty. Residue indices are sorted and unique across both chains.
///
/// Runs in `O(chain atoms · local density)` time.
///
/// # Errors
///
/// Returns [`InterfaceError`] for a non-finite or negative cutoff, an ambiguous
/// name, or two selectors that resolve to the same chain.
pub fn chain_interface(
    structure: &Structure,
    first: &str,
    second: &str,
    cutoff: f32,
    backend: SpatialBackend,
    context: &ExecutionContext,
) -> Result<Vec<ResidueIndex>, InterfaceError> {
    chain_interface_with_options(
        structure,
        &InterfaceOptions {
            first: ChainSelector::any(first),
            second: ChainSelector::any(second),
            cutoff,
            backend,
            periodic: false,
        },
        context,
    )
}

/// Returns the interface residues for explicit chain selectors.
///
/// With `periodic` set the structure's unit cell supplies the minimum-image
/// convention, so contacts across a cell boundary count.
///
/// # Errors
///
/// Returns [`InterfaceError`] for a bad cutoff, an ambiguous or identical pair
/// of selectors, or a missing or placeholder cell when `periodic` is set.
pub fn chain_interface_with_options(
    structure: &Structure,
    options: &InterfaceOptions<'_>,
    context: &ExecutionContext,
) -> Result<Vec<ResidueIndex>, InterfaceError> {
    let (first_atoms, second_atoms) = resolve_pair(structure, options)?;
    let periodic_box = periodic_box(structure, options.periodic)?;
    let left = AtomSelection::from_sorted(first_atoms);
    let right = AtomSelection::from_sorted(second_atoms);
    let query = PairQuery {
        positions: structure.positions(),
        left: &left,
        right: &right,
        cutoff: options.cutoff,
        options: SpatialSearchOptions::with_backend(options.backend),
        periodic: periodic_box.as_ref(),
        context,
    };
    let parts = reduce_pairs_within_unsorted(&query, Vec::new, |residues: &mut Vec<u32>, pair| {
        push_interface_residues(structure, pair, residues);
    })?;

    let mut residues: Vec<u32> = Vec::new();
    for part in parts {
        residues.extend(part);
    }

    Ok(finish_interface_residues(residues))
}

/// Returns interface residues through a structure-bound spatial resolver.
///
/// A compiled plan uses this entry point so repeated interface/contact
/// operations can share the resolver's bounded cell/k-d index cache.  The
/// resolver remains borrowed from the same immutable structure snapshot, which
/// makes the lifetime and coordinate generation relationship explicit.
///
/// # Errors
///
/// Returns a diagnostic when the resolver rejects the workload, a chain name is
/// ambiguous or both names resolve to the same chain.
pub fn chain_interface_with_spatial(
    structure: &Structure,
    first: &str,
    second: &str,
    cutoff: f32,
    backend: SpatialBackend,
    spatial: &StructureSpatial<'_>,
) -> Result<Vec<ResidueIndex>, Diagnostic> {
    let options = InterfaceOptions {
        first: ChainSelector::any(first),
        second: ChainSelector::any(second),
        cutoff,
        backend,
        periodic: false,
    };
    let (first_atoms, second_atoms) = resolve_pair(structure, &options)?;
    let left = AtomSelection::from_sorted(first_atoms);
    let right = AtomSelection::from_sorted(second_atoms);
    let pairs = spatial.pairs_with_backend(&left, &right, cutoff, backend)?;

    Ok(interface_residues(structure, pairs))
}

fn periodic_box(
    structure: &Structure,
    periodic: bool,
) -> Result<Option<PeriodicBox>, InterfaceError> {
    if !periodic {
        return Ok(None);
    }
    let cell = structure.data().cell.ok_or(InterfaceError::MissingCell)?;
    if cell.is_placeholder() {
        return Err(InterfaceError::PlaceholderCell);
    }
    Ok(Some(PeriodicBox::from_cell(cell)?))
}

fn interface_residues(
    structure: &Structure,
    pairs: Vec<molframe_spatial::NeighborPair>,
) -> Vec<ResidueIndex> {
    let mut residues = Vec::new();
    for pair in pairs {
        push_interface_residues(structure, pair, &mut residues);
    }
    finish_interface_residues(residues)
}

/// Records both endpoints' residues for one contact pair.
///
/// Splitting the accumulation from the finish lets a streaming query feed the
/// same reduction without first materialising a pair vector.
fn push_interface_residues(
    structure: &Structure,
    pair: molframe_spatial::NeighborPair,
    residues: &mut Vec<u32>,
) {
    if let Some(residue) = residue_of(structure, pair.first) {
        residues.push(residue);
    }
    if let Some(residue) = residue_of(structure, pair.second) {
        residues.push(residue);
    }
}

/// Orders and deduplicates accumulated residues.
///
/// The sort makes the result independent of the order pairs arrived in, so a
/// streaming query and a materialising one agree exactly.
fn finish_interface_residues(mut residues: Vec<u32>) -> Vec<ResidueIndex> {
    residues.sort_unstable();
    residues.dedup();
    residues.into_iter().map(ResidueIndex::new).collect()
}

/// The chains a selector denotes, as chain indices.
fn matching_chains(
    structure: &Structure,
    selector: ChainSelector<'_>,
) -> Result<Vec<u32>, InterfaceError> {
    let by_label: Vec<u32> = chain_indices(structure, |chain| chain.label() == Some(selector.name));
    let by_auth: Vec<u32> =
        chain_indices(structure, |chain| chain.auth_label() == Some(selector.name));
    match selector.namespace {
        Some(ChainNamespace::Label) => Ok(by_label),
        Some(ChainNamespace::Auth) => Ok(by_auth),
        None if by_label.is_empty() => Ok(by_auth),
        None if by_auth.is_empty() || by_label == by_auth => Ok(by_label),
        None => Err(InterfaceError::AmbiguousChain(selector.name.to_owned())),
    }
}

fn chain_indices(structure: &Structure, keep: impl Fn(ChainRef<'_>) -> bool) -> Vec<u32> {
    structure
        .data()
        .chains()
        .filter(|chain| keep(*chain))
        .map(|chain| chain.index().get())
        .collect()
}

/// Resolves both selectors to atom lists, rejecting a shared chain.
fn resolve_pair(
    structure: &Structure,
    options: &InterfaceOptions<'_>,
) -> Result<(Vec<u32>, Vec<u32>), InterfaceError> {
    let first = matching_chains(structure, options.first)?;
    let second = matching_chains(structure, options.second)?;
    if first.iter().any(|chain| second.contains(chain)) {
        return Err(InterfaceError::SameChain);
    }
    Ok((
        chain_atoms(structure, &first),
        chain_atoms(structure, &second),
    ))
}

/// Collects the sorted atom indices of the given chains.
fn chain_atoms(structure: &Structure, chains: &[u32]) -> Vec<u32> {
    let mut atoms = Vec::new();
    for chain in structure.data().chains() {
        if !chains.contains(&chain.index().get()) {
            continue;
        }
        for residue in chain.residues() {
            for atom in residue.atoms() {
                atoms.push(atom.index().get());
            }
        }
    }
    atoms.sort_unstable();
    atoms
}

fn residue_of(structure: &Structure, atom: u32) -> Option<u32> {
    structure
        .data()
        .atom(molframe_core::index::AtomIndex::new(atom))
        .and_then(AtomRef::residue)
        .map(|residue| residue.index().get())
}

#[cfg(test)]
#[path = "interface_tests.rs"]
mod tests;
