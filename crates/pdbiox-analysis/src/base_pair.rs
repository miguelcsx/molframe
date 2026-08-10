//! Watson–Crick pairing from CCD base identity and oriented hydrogen bonds.

use crate::{HydrogenBondError, HydrogenBondOptions, hydrogen_bonds};
use pdbiox_chem::{ComponentKind, ComponentProvider};
use pdbiox_core::index::ResidueIndex;
use pdbiox_core::{AtomAnnotation, Diagnostic, Presence, Structure};
use std::collections::BTreeMap;

/// Explicit chemical and geometric policy for canonical base pairing.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BasePairOptions {
    /// Hydrogen-bond geometry and periodic policy.
    pub hydrogen_bonds: HydrogenBondOptions,
    /// Minimum number of oriented inter-base hydrogen bonds.
    pub minimum_hydrogen_bonds: usize,
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

/// Finds canonical pairs supported by CCD identity and oriented hydrogen bonds.
///
/// Modified nucleotides participate when their CCD entry supplies a canonical
/// one-letter code. No residue names or atom-name edge tables are embedded.
///
/// # Errors
///
/// Returns [`BasePairError`] for invalid policy, missing CCD nucleotide data,
/// provider failure, or unavailable hydrogen-bond chemistry.
pub fn base_pairs(
    structure: &Structure,
    provider: &dyn ComponentProvider,
    options: BasePairOptions,
) -> Result<Vec<BasePair>, BasePairError> {
    if options.minimum_hydrogen_bonds == 0 {
        return Err(BasePairError::InvalidOptions);
    }
    let bases = base_identities(structure, provider)?;
    let bonds = hydrogen_bonds(structure, options.hydrogen_bonds)?;
    let mut support: BTreeMap<(ResidueIndex, ResidueIndex), Vec<f32>> = BTreeMap::new();
    for bond in bonds {
        let Some(donor_residue) = structure
            .atom(bond.donor)
            .and_then(pdbiox_core::structure::AtomRef::residue)
        else {
            continue;
        };
        let Some(acceptor_residue) = structure
            .atom(bond.acceptor)
            .and_then(pdbiox_core::structure::AtomRef::residue)
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
        let pair = ordered(donor_residue.index(), acceptor_residue.index());
        support
            .entry(pair)
            .or_default()
            .push(bond.donor_acceptor_distance);
    }
    Ok(support
        .into_iter()
        .filter(|(_, distances)| distances.len() >= options.minimum_hydrogen_bonds)
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
        .get(pdbiox_core::COMPONENT_KIND_ANNOTATION)?
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
