//! Explicit distance-based connectivity inference.

use crate::{
    AtomIndex, AtomSelection, BondOrder, BondProvenance, BondRecord, Code, Diagnostic,
    SpatialBackend, Structure,
};
use pdbiox_core::bond::BondTableBuilder;
use std::collections::BTreeSet;

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
            scale: 1.15,
            lower_bound: 0.4,
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
) -> Result<BondInferenceReport, Diagnostic> {
    if !options.scale.is_finite()
        || options.scale <= 0.0
        || !options.lower_bound.is_finite()
        || options.lower_bound < 0.0
    {
        return Err(Diagnostic::new(Code::E4002).with_context("parameter", "bond inference"));
    }
    if !structure.data().coords.is_dense() {
        return Err(Diagnostic::new(Code::E6008));
    }
    let radii: Vec<_> = structure
        .data()
        .atoms()
        .map(|atom| {
            atom.element()
                .and_then(pdbiox_chem::element_properties)
                .and_then(|properties| properties.covalent_radius)
        })
        .collect();
    let skipped_atoms = radii
        .iter()
        .enumerate()
        .filter(|(_, radius)| radius.is_none())
        .map(|(atom, _)| AtomIndex::new(atom as u32))
        .collect();
    let max_radius = radii.iter().flatten().copied().fold(0.0f32, f32::max);
    let maximum = max_radius * 2.0 * options.scale;
    let pairs = pdbiox_spatial::pairs_within(
        structure.positions(),
        &AtomSelection::All(structure.atom_count()),
        &AtomSelection::All(structure.atom_count()),
        maximum,
        options.backend,
        None,
    )
    .map_err(|error| Diagnostic::new(Code::E9001).with_message(error.to_string()))?;
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
    for pair in pairs {
        if reject_pair(&pair, &radii, &chains, &existing, options) {
            continue;
        }
        output.push(BondRecord {
            atom_a: AtomIndex::new(pair.first),
            atom_b: AtomIndex::new(pair.second),
            order: BondOrder::Single,
            provenance: BondProvenance::InferredDistance,
        });
    }
    let mut data = structure.data().clone();
    data.bonds = output.finish();
    Ok(BondInferenceReport {
        structure: Structure::new(data),
        skipped_atoms,
    })
}

fn reject_pair(
    pair: &pdbiox_spatial::NeighborPair,
    radii: &[Option<f32>],
    chains: &[u32],
    existing: &BTreeSet<(u32, u32)>,
    options: BondInference,
) -> bool {
    if options.respect_existing && existing.contains(&(pair.first, pair.second)) {
        return true;
    }
    if options.exclude_across_chains
        && chains.get(pair.first as usize) != chains.get(pair.second as usize)
    {
        return true;
    }
    let (Some(Some(left)), Some(Some(right))) = (
        radii.get(pair.first as usize),
        radii.get(pair.second as usize),
    ) else {
        return true;
    };
    let lower = options.lower_bound * options.lower_bound;
    let upper = (left + right).powi(2) * options.scale.powi(2);
    pair.distance_squared < lower || pair.distance_squared > upper
}

fn atom_chains(structure: &Structure) -> Vec<u32> {
    let mut chains = vec![u32::MAX; structure.atom_count() as usize];
    for chain in structure.data().chains() {
        for residue in chain.residues() {
            for atom in residue.atoms() {
                chains[atom.index().as_usize()] = chain.index().get();
            }
        }
    }
    chains
}

#[cfg(test)]
#[path = "chemistry_api_tests.rs"]
mod tests;
