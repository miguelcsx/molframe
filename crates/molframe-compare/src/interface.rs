//! Chain atoms and inter-set contacts, shared by the interface-based scores.
//!
//! Both `DockQ` and the quaternary-structure score start from the same two
//! questions: which atoms belong to a chain, and which pairs across two atom
//! sets are in contact. Answering them once here keeps the two scores agreeing
//! on what an interface contact is. The residue-level scores also share the
//! per-atom annotations built here: which residue an atom is in, whether it is a
//! heavy atom, whether it is a backbone atom, and whether it stands for its
//! residue in side-chain-centroid style contact maps.

use molframe_core::contract::Namespace;
use molframe_core::element::Element;
use molframe_core::structure::Structure;

use crate::CompareError;

/// The atom indices belonging to a chain in one explicit namespace, sorted.
pub(crate) fn chain_atoms(
    structure: &Structure,
    name: &str,
    namespace: Namespace,
) -> Result<Vec<usize>, CompareError> {
    let mut atoms = Vec::new();
    for chain in structure.data().chains() {
        let matches = match namespace {
            Namespace::Label => chain.label() == Some(name),
            Namespace::Auth => chain.auth_label() == Some(name),
            _ => return Err(CompareError::UnsupportedNamespace),
        };
        if !matches {
            continue;
        }
        for residue in chain.residues() {
            for atom in residue.atoms() {
                atoms.push(atom.index().as_usize());
            }
        }
    }
    atoms.sort_unstable();
    Ok(atoms)
}

/// The pairs `(a, b)` from the two index sets whose atoms are within `cutoff`.
pub(crate) fn contacts(
    positions: &[[f32; 3]],
    first: &[usize],
    second: &[usize],
    cutoff: f32,
) -> Vec<(usize, usize)> {
    let cutoff_squared = f64::from(cutoff) * f64::from(cutoff);
    let mut pairs = Vec::new();
    for &a in first {
        let Some(&point_a) = positions.get(a) else {
            continue;
        };
        for &b in second {
            let Some(&point_b) = positions.get(b) else {
                continue;
            };
            if molframe_geom::distance_squared(point_a, point_b) <= cutoff_squared {
                pairs.push((a, b));
            }
        }
    }
    pairs
}

/// What a residue-level score needs to know about one atom.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct AtomInfo {
    /// Structure-wide residue number; unique across chains.
    pub residue: usize,
    /// Not a hydrogen or deuterium.
    pub heavy: bool,
    /// One of the protein backbone atoms N, CA, C, O.
    pub backbone: bool,
    /// The residue's representative atom: CB, or CA for glycine.
    pub representative: bool,
}

/// Annotation for an atom that belongs to no residue; it takes part in nothing.
pub(crate) const UNANNOTATED: AtomInfo = AtomInfo {
    residue: usize::MAX,
    heavy: false,
    backbone: false,
    representative: false,
};

/// The annotation of every atom of a structure, indexed by atom index.
pub(crate) fn atom_infos(structure: &Structure) -> Vec<AtomInfo> {
    let mut infos = vec![UNANNOTATED; structure.atom_count() as usize];
    let mut residue_number = 0usize;
    for chain in structure.data().chains() {
        for residue in chain.residues() {
            let glycine = residue.name() == Some("GLY");
            let mut beta: Option<usize> = None;
            let mut alpha: Option<usize> = None;
            for atom in residue.atoms() {
                let index = atom.index().as_usize();
                let name = match atom.name() {
                    Some(name) => name,
                    None => "",
                };
                let element = atom.element();
                if name == "CB" && beta.is_none() && element_fits(element, Element::CARBON) {
                    beta = Some(index);
                }
                if name == "CA" && alpha.is_none() && element_fits(element, Element::CARBON) {
                    alpha = Some(index);
                }
                if let Some(slot) = infos.get_mut(index) {
                    *slot = AtomInfo {
                        residue: residue_number,
                        heavy: is_heavy(name, element),
                        backbone: is_backbone(name, element),
                        representative: false,
                    };
                }
            }
            let representative = match (beta, alpha) {
                (Some(index), _) => Some(index),
                (None, Some(index)) if glycine => Some(index),
                _ => None,
            };
            if let Some(slot) = representative.and_then(|index| infos.get_mut(index)) {
                slot.representative = true;
            }
            residue_number += 1;
        }
    }
    infos
}

/// An unknown element fits anything; a known one must match.
fn element_fits(element: Option<Element>, expected: Element) -> bool {
    match element {
        Some(found) if !found.is_unknown() => found == expected,
        _ => true,
    }
}

fn is_heavy(name: &str, element: Option<Element>) -> bool {
    let element = match element {
        Some(found) if !found.is_unknown() => found,
        _ => Element::infer_from_name(name),
    };
    !element.is_hydrogen()
}

fn is_backbone(name: &str, element: Option<Element>) -> bool {
    match name {
        "N" => element_fits(element, Element::NITROGEN),
        "CA" | "C" => element_fits(element, Element::CARBON),
        "O" => element_fits(element, Element::OXYGEN),
        _ => false,
    }
}
