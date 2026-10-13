//! Which atoms get a radius, and from where, for surface and steric work.
//!
//! A structure file mixes the molecule of interest with crystallisation water,
//! ligands and ions. Whether those belong in a surface or contact calculation
//! is a question the caller must be able to answer explicitly, and the default
//! must be the one the field's reference programs share: NACCESS and `FreeSASA`
//! both measure the polymer alone and ignore heteroatoms unless told otherwise.
//! An excluded atom gets no radius (`NaN`), so it neither has an area nor
//! occludes its neighbours.
//!
//! Metals are handled by one rule. A radius set that has no value for a metal
//! (Bondi has none for Ca, Zn, Fe, Cu, Ni, Co or Mn) either leaves the atom
//! without a radius or draws on the metal table, as the caller chooses; the
//! choice is never made silently.
//!
//! One pass over the residues, `O(N)` in atoms, one allocation.

use crate::metal::{MetalRadiusKind, metal_radius};
use crate::standard_bonds::{is_amino_acid_component, is_water_component};
use crate::{ComponentKind, RadiusSet, vdw_radius};
use molframe_core::annotation::AtomAnnotation;
use molframe_core::column::Presence;
use molframe_core::structure::{ResidueRef, Structure};

/// What a residue is, for the purpose of deciding whether it is measured.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ResidueClass {
    /// Polymer residue, including modified amino acids written as heterogens.
    Polymer,
    /// Water.
    Water,
    /// A monatomic ion.
    Ion,
    /// Any other non-polymer: cofactors, drugs, buffers, sugars.
    Ligand,
}

/// Which heteroatom classes take part. The default measures the polymer only.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct HeteroPolicy {
    /// Include water.
    pub waters: bool,
    /// Include ligands.
    pub ligands: bool,
    /// Include monatomic ions.
    pub ions: bool,
}

impl HeteroPolicy {
    /// Polymer only: the reference-program default.
    pub const POLYMER: Self = Self {
        waters: false,
        ligands: false,
        ions: false,
    };
    /// Polymer, ligands and ions, without water.
    pub const SOLUTE: Self = Self {
        waters: false,
        ligands: true,
        ions: true,
    };
    /// Everything in the file.
    pub const ALL: Self = Self {
        waters: true,
        ligands: true,
        ions: true,
    };

    /// Whether residues of `class` take part.
    #[must_use]
    pub const fn includes(self, class: ResidueClass) -> bool {
        match class {
            ResidueClass::Polymer => true,
            ResidueClass::Water => self.waters,
            ResidueClass::Ion => self.ions,
            ResidueClass::Ligand => self.ligands,
        }
    }
}

/// What to do for a metal the chosen radius set has no value for.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum MetalFallback {
    /// Leave the atom without a radius.
    None,
    /// Use the Alvarez van der Waals radius: the ion as a solvent obstacle.
    #[default]
    Vdw,
    /// Use the Shannon ionic radius: the ion as a steric contact.
    Ionic,
}

/// A complete, explicit radius assignment.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct RadiusPolicy {
    /// Radius set for every element it covers.
    pub set: RadiusSet,
    /// Source for metals the set does not cover.
    pub metals: MetalFallback,
    /// Heteroatom classes that take part.
    pub hetero: HeteroPolicy,
}

impl Default for RadiusPolicy {
    fn default() -> Self {
        Self {
            set: RadiusSet::Bondi,
            metals: MetalFallback::default(),
            hetero: HeteroPolicy::POLYMER,
        }
    }
}

/// The radius of one element from `set`, drawing on the metal table under
/// `metals` when the set has no value for it.
#[must_use]
pub fn radius_with_fallback(
    element: molframe_core::Element,
    set: RadiusSet,
    metals: MetalFallback,
) -> Option<f32> {
    vdw_radius(element, set).or_else(|| match metals {
        MetalFallback::None => None,
        MetalFallback::Vdw => metal_radius(element, MetalRadiusKind::Vdw),
        MetalFallback::Ionic => metal_radius(element, MetalRadiusKind::Ionic),
    })
}

/// The class of one residue.
///
/// A component-kind column, when present, decides; otherwise the residue name
/// and the heterogen flag do. A heterogen with a single atom is an ion.
#[must_use]
pub fn classify_residue(residue: &ResidueRef<'_>, kind: Option<ComponentKind>) -> ResidueClass {
    let name = residue.name();
    if name.is_some_and(is_water_component) {
        return ResidueClass::Water;
    }
    match kind {
        Some(ComponentKind::AminoAcid | ComponentKind::Nucleotide) => {
            return ResidueClass::Polymer;
        }
        Some(ComponentKind::Solvent) => return ResidueClass::Water,
        Some(ComponentKind::Ion) => return ResidueClass::Ion,
        Some(ComponentKind::Saccharide | ComponentKind::Lipid | ComponentKind::NonPolymer) => {
            return ResidueClass::Ligand;
        }
        Some(ComponentKind::Unknown) | None => {}
    }
    if !residue.is_het() || name.is_some_and(is_amino_acid_component) {
        ResidueClass::Polymer
    } else if residue.atoms().take(2).count() == 1 {
        ResidueClass::Ion
    } else {
        ResidueClass::Ligand
    }
}

/// The radius of every atom of `structure` under `policy`, indexed by atom.
///
/// `NaN` marks an atom that is excluded by the heteroatom policy, whose element
/// is unknown, or which has no value in the set and no metal fallback.
#[must_use]
pub fn policy_radii(structure: &Structure, policy: &RadiusPolicy) -> Vec<f32> {
    let mut radii = vec![f32::NAN; structure.atom_count() as usize];
    let kinds = match structure
        .annotations()
        .get(molframe_core::COMPONENT_KIND_ANNOTATION)
    {
        Some(AtomAnnotation::Integer(column)) => Some(column),
        _ => None,
    };
    for chain in structure.data().chains() {
        for residue in chain.residues() {
            let declared = kinds.and_then(|column| {
                let first = residue.atoms().next()?;
                match column.get(first.index().get()) {
                    Some((code, Presence::Present)) => ComponentKind::from_code(code),
                    _ => None,
                }
            });
            if !policy.hetero.includes(classify_residue(&residue, declared)) {
                continue;
            }
            for atom in residue.atoms() {
                let radius = atom
                    .element()
                    .and_then(|element| radius_with_fallback(element, policy.set, policy.metals));
                if let (Some(radius), Some(slot)) = (radius, radii.get_mut(atom.index().as_usize()))
                {
                    *slot = radius;
                }
            }
        }
    }
    radii
}

#[cfg(test)]
#[path = "radius_policy_tests.rs"]
mod tests;
