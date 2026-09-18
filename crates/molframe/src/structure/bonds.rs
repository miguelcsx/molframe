//! Explicit distance-based connectivity inference.

use crate::{
    AtomIndex, AtomSelection, BondOrder, BondProvenance, BondRecord, Code, Diagnostic, Findings,
    SpatialBackend, Structure,
};
use molframe_core::bond::BondTableBuilder;
use std::collections::BTreeSet;

/// Default multiplier applied to the sum of two CCD covalent radii.
pub const DEFAULT_BOND_RADIUS_SCALE: f32 = 1.15;

/// Default minimum accepted interatomic distance in angstroms.
pub const DEFAULT_MINIMUM_BOND_DISTANCE: f32 = 0.4;

const BOND_INFERENCE_PARAMETER: &str = "bond inference";

/// Parameters for explicitly requested geometric bond inference.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BondInference {
    /// Multiplier on the sum of covalent radii.
    pub scale: f32,
    /// Distances below this bound are rejected as atomic overlap.
    pub lower_bound: f32,
    /// Whether candidates may cross chain boundaries.
    pub exclude_across_chains: bool,
    /// Whether existing connectivity has priority over inferred candidates.
    pub respect_existing: bool,
    /// Spatial implementation, with `Auto` selecting by workload.
    pub backend: SpatialBackend,
}

impl Default for BondInference {
    fn default() -> Self {
        Self {
            scale: DEFAULT_BOND_RADIUS_SCALE,
            lower_bound: DEFAULT_MINIMUM_BOND_DISTANCE,
            exclude_across_chains: false,
            respect_existing: true,
            backend: SpatialBackend::Auto,
        }
    }
}

/// A new connectivity snapshot and atoms lacking a covalent radius.
#[derive(Clone, Debug)]
pub struct BondInferenceReport {
    /// Structure with inferred edges added.
    pub structure: Structure,
    /// Atoms excluded because their element has no radius in the table.
    pub skipped_atoms: Vec<AtomIndex>,
}

/// Infers bonds by covalent-radius distance using the shared spatial planner.
///
/// The operation is never implicit during parsing. Existing records are copied
/// first, so file and CCD provenance cannot be overwritten by an inference.
///
/// # Errors
///
/// Returns a diagnostic for invalid bounds, ragged coordinates, or a spatial
/// indexing failure.
pub fn infer_bonds(
    structure: &Structure,
    options: BondInference,
    context: &molframe_core::ExecutionContext,
) -> Result<BondInferenceReport, Findings> {
    if !options.scale.is_finite()
        || options.scale <= 0.0
        || !options.lower_bound.is_finite()
        || options.lower_bound < 0.0
    {
        return Err(Diagnostic::new(Code::E4002)
            .with_context("parameter", BOND_INFERENCE_PARAMETER)
            .into());
    }
    if !structure.data().coords.is_dense() {
        return Err(Diagnostic::new(Code::E6008).into());
    }
    let radii: Vec<_> = structure
        .data()
        .atoms()
        .map(|atom| {
            atom.element()
                .and_then(molframe_chem::element_properties)
                .and_then(|properties| properties.covalent_radius)
        })
        .collect();
    let skipped_atoms = radii
        .iter()
        .enumerate()
        .filter(|(_, radius)| radius.is_none())
        .map(|(atom, _)| {
            u32::try_from(atom)
                .map(AtomIndex::new)
                .map_err(|error| Diagnostic::new(Code::E9001).with_message(error.to_string()))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let max_radius = radii.iter().flatten().copied().fold(0.0f32, f32::max);
    let maximum = max_radius * 2.0 * options.scale;
    let lower_squared = options.lower_bound * options.lower_bound;
    if !maximum.is_finite() || !lower_squared.is_finite() || options.lower_bound > maximum {
        return Err(Diagnostic::new(Code::E4002)
            .with_context("parameter", BOND_INFERENCE_PARAMETER)
            .into());
    }
    let chains = atom_chains(structure);
    let existing: BTreeSet<_> = structure
        .data()
        .bonds
        .iter()
        .map(|bond| (bond.atom_a.get(), bond.atom_b.get()))
        .collect();
    let mut output = BondTableBuilder::new();
    for bond in structure.data().bonds.iter() {
        output.push(bond);
    }
    // Only pairs that become bonds survive, a small fraction of the candidates,
    // so the candidates are filtered as they are produced rather than collected
    // into a vector sized by the quadratic candidate count.
    let selection = AtomSelection::All(structure.atom_count());
    molframe_spatial::for_each_pairs_within_unsorted(
        &molframe_spatial::PairQuery {
            positions: structure.positions(),
            left: &selection,
            right: &selection,
            cutoff: maximum,
            options: molframe_spatial::SpatialSearchOptions::with_backend(options.backend),
            periodic: None,
            context,
        },
        |pair| {
            if reject_pair(&pair, &radii, &chains, &existing, options, lower_squared) {
                return;
            }
            output.push(BondRecord {
                atom_a: AtomIndex::new(pair.first),
                atom_b: AtomIndex::new(pair.second),
                order: BondOrder::Single,
                provenance: BondProvenance::InferredDistance,
            });
        },
    )
    .map_err(|error| Diagnostic::new(Code::E9001).with_message(error.to_string()))?;

    let mut data = structure.data().clone();
    data.bonds = output.finish();
    Ok(BondInferenceReport {
        structure: Structure::new(data),
        skipped_atoms,
    })
}

fn reject_pair(
    pair: &molframe_spatial::NeighborPair,
    radii: &[Option<f32>],
    chains: &[Option<u32>],
    existing: &BTreeSet<(u32, u32)>,
    options: BondInference,
    lower_squared: f32,
) -> bool {
    if options.respect_existing && existing.contains(&(pair.first, pair.second)) {
        return true;
    }
    if options.exclude_across_chains
        && chains.get(AtomIndex::new(pair.first).as_usize())
            != chains.get(AtomIndex::new(pair.second).as_usize())
    {
        return true;
    }
    let (Some(Some(left)), Some(Some(right))) = (
        radii.get(AtomIndex::new(pair.first).as_usize()),
        radii.get(AtomIndex::new(pair.second).as_usize()),
    ) else {
        return true;
    };
    let upper = (left + right).powi(2) * options.scale.powi(2);
    pair.distance_squared < lower_squared || pair.distance_squared > upper
}

fn atom_chains(structure: &Structure) -> Vec<Option<u32>> {
    let mut chains = vec![None; AtomIndex::new(structure.atom_count()).as_usize()];
    for chain in structure.data().chains() {
        for residue in chain.residues() {
            for atom in residue.atoms() {
                chains[atom.index().as_usize()] = Some(chain.index().get());
            }
        }
    }
    chains
}

#[cfg(test)]
#[path = "bonds_tests.rs"]
mod tests;
