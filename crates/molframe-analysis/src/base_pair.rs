//! Watson–Crick pairing from CCD base identity, oriented hydrogen bonds between
//! the Watson–Crick edge atoms, and base reference geometry.

use crate::{HydrogenBondError, HydrogenBondOptions, hydrogen_bonds};
use molframe_chem::{ComponentKind, ComponentProvider};
use molframe_core::index::ResidueIndex;
use molframe_core::{AtomAnnotation, Diagnostic, ExecutionContext, Presence, Structure};
use std::collections::BTreeMap;

#[path = "base_pair_edges.rs"]
mod edges;

pub use edges::WatsonCrickGeometry;

/// Explicit chemical and geometric policy for canonical base pairing.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BasePairOptions {
    /// Hydrogen-bond geometry and periodic policy.
    pub hydrogen_bonds: HydrogenBondOptions,
    /// Minimum number of oriented inter-base hydrogen bonds.
    pub minimum_hydrogen_bonds: usize,
    /// Reference-geometry limits (C1'–C1' distance and base-plane alignment)
    /// that a pair must meet to count as Watson–Crick. `None` skips the
    /// geometry test and keeps only the edge-atom test; hydrogen bonds that do
    /// not join Watson–Crick edge atoms never count either way.
    pub geometry: Option<WatsonCrickGeometry>,
}

/// A Watson–Crick base pair between two CCD-identified nucleotide residues.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BasePair {
    /// Lower-indexed residue.
    pub first: ResidueIndex,
    /// Higher-indexed residue.
    pub second: ResidueIndex,
    /// Number of oriented hydrogen bonds supporting the pair.
    pub hydrogen_bond_count: usize,
    /// Shortest donor–acceptor distance among the supporting bonds.
    pub closest_distance: f32,
}

/// Why base-pair detection could not be evaluated.
#[derive(Debug, thiserror::Error)]
pub enum BasePairError {
    /// At least one supporting hydrogen bond must be requested.
    #[error("minimum hydrogen-bond count must be non-zero")]
    InvalidOptions,
    /// A nucleotide annotation names a component absent from the explicit CCD.
    #[error("CCD component {component} for residue {residue} is unavailable")]
    MissingComponent {
        /// Residue whose component could not be resolved.
        residue: ResidueIndex,
        /// Deposited component identifier.
        component: Box<str>,
    },
    /// Component provider access failed.
    #[error("component provider failed: {0}")]
    Provider(Diagnostic),
    /// Hydrogen-bond chemistry or geometry failed.
    #[error(transparent)]
    HydrogenBond(#[from] HydrogenBondError),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum CanonicalBase {
    Adenine,
    Cytosine,
    Guanine,
    Thymine,
    Uracil,
}

impl CanonicalBase {
    const fn from_one_letter(code: u8) -> Option<Self> {
        match code.to_ascii_uppercase() {
            b'A' => Some(Self::Adenine),
            b'C' => Some(Self::Cytosine),
            b'G' => Some(Self::Guanine),
            b'T' => Some(Self::Thymine),
            b'U' => Some(Self::Uracil),
            _ => None,
        }
    }

    const fn complementary(self, other: Self) -> bool {
        matches!(
            (self, other),
            (Self::Adenine, Self::Thymine | Self::Uracil)
                | (Self::Thymine | Self::Uracil, Self::Adenine)
                | (Self::Guanine, Self::Cytosine)
                | (Self::Cytosine, Self::Guanine)
        )
    }
}

/// Finds canonical pairs supported by CCD identity, Watson–Crick-edge hydrogen
/// bonds and, when requested, reference geometry.
///
/// Modified nucleotides participate when their CCD entry supplies a canonical
/// one-letter code. Only hydrogen bonds joining the Watson–Crick edge atoms
/// (A–T/U: N6…O4, N1…N3; G–C: O6…N4, N1…N3, N2…O2) support a pair, so a
/// Hoogsteen or sugar-edge contact between complementary bases is left out.
/// Candidate pairs must also meet the [`WatsonCrickGeometry`] limits.
///
/// # Errors
///
/// Returns [`BasePairError`] for invalid policy, missing CCD nucleotide data,
/// provider failure, or unavailable hydrogen-bond chemistry.
pub fn base_pairs(
    structure: &Structure,
    provider: &dyn ComponentProvider,
    options: BasePairOptions,
    context: &ExecutionContext,
) -> Result<Vec<BasePair>, BasePairError> {
    if options.minimum_hydrogen_bonds == 0 {
        return Err(BasePairError::InvalidOptions);
    }
    let bases = base_identities(structure, provider)?;
    let bonds = hydrogen_bonds(structure, options.hydrogen_bonds, context)?;
    let mut support: BTreeMap<(ResidueIndex, ResidueIndex), Vec<f32>> = BTreeMap::new();
    for bond in bonds.iter() {
        let Some(donor_residue) = structure
            .atom(bond.donor)
            .and_then(molframe_core::structure::AtomRef::residue)
        else {
            continue;
        };
        let Some(acceptor_residue) = structure
            .atom(bond.acceptor)
            .and_then(molframe_core::structure::AtomRef::residue)
        else {
            continue;
        };
        if donor_residue.index() == acceptor_residue.index() {
            continue;
        }
        let Some(&donor_base) = bases.get(&donor_residue.index()) else {
            continue;
        };
        let Some(&acceptor_base) = bases.get(&acceptor_residue.index()) else {
            continue;
        };
        if !donor_base.complementary(acceptor_base) {
            continue;
        }
        let (Some(donor_name), Some(acceptor_name)) = (
            structure
                .atom(bond.donor)
                .and_then(molframe_core::structure::AtomRef::name),
            structure
                .atom(bond.acceptor)
                .and_then(molframe_core::structure::AtomRef::name),
        ) else {
            continue;
        };
        if !edges::watson_crick_contact(donor_base, donor_name, acceptor_base, acceptor_name) {
            continue;
        }
        let pair = ordered(donor_residue.index(), acceptor_residue.index());
        support
            .entry(pair)
            .or_default()
            .push(bond.donor_acceptor_distance);
    }
    Ok(support
        .into_iter()
        .filter(|(_, distances)| distances.len() >= options.minimum_hydrogen_bonds)
        .filter(|((first, second), _)| match options.geometry {
            None => true,
            Some(limits) => {
                match (
                    structure.data().residue(*first),
                    structure.data().residue(*second),
                ) {
                    (Some(a), Some(b)) => edges::has_watson_crick_geometry(a, b, limits),
                    _ => false,
                }
            }
        })
        .map(|((first, second), distances)| BasePair {
            first,
            second,
            hydrogen_bond_count: distances.len(),
            closest_distance: distances.into_iter().fold(f32::INFINITY, f32::min),
        })
        .collect())
}

fn base_identities(
    structure: &Structure,
    provider: &dyn ComponentProvider,
) -> Result<BTreeMap<ResidueIndex, CanonicalBase>, BasePairError> {
    let mut output = BTreeMap::new();
    for residue in structure.data().residues() {
        if residue_component_kind(structure, residue.index()) != Some(ComponentKind::Nucleotide) {
            continue;
        }
        let Some(component_id) = residue.name() else {
            continue;
        };
        let component = provider
            .get(component_id)
            .map_err(BasePairError::Provider)?
            .ok_or_else(|| BasePairError::MissingComponent {
                residue: residue.index(),
                component: component_id.into(),
            })?;
        if let Some(base) = component
            .one_letter_code
            .and_then(CanonicalBase::from_one_letter)
        {
            output.insert(residue.index(), base);
        }
    }
    Ok(output)
}

fn residue_component_kind(structure: &Structure, residue: ResidueIndex) -> Option<ComponentKind> {
    let AtomAnnotation::Integer(column) = structure
        .annotations()
        .get(molframe_core::COMPONENT_KIND_ANNOTATION)?
    else {
        return None;
    };
    structure.data().residue(residue)?.atoms().find_map(|atom| {
        column
            .get(atom.index().get())
            .filter(|(_, presence)| *presence == Presence::Present)
            .and_then(|(code, _)| ComponentKind::from_code(code))
    })
}

fn ordered(first: ResidueIndex, second: ResidueIndex) -> (ResidueIndex, ResidueIndex) {
    if first <= second {
        (first, second)
    } else {
        (second, first)
    }
}

#[cfg(test)]
#[path = "base_pair_tests.rs"]
mod tests;
