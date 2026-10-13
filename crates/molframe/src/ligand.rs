//! Comparing poses of one ligand held as structures: the heavy-atom graph is read off the
//! structure and handed to the symmetry-aware RMSD of [`crate::compare`], so a flipped phenyl
//! ring is not a different pose.

use crate::structure::Structure;
use crate::{BondOrder, compare::ligand_symmetry_rmsd};
use molframe_chem::{Component, ComponentAtom, ComponentBond, ComponentKind};
use std::collections::HashMap;
use std::sync::Arc;

/// The heavy atoms of a ligand as a chemical component, with the position in the ligand
/// structure that each component atom stands for.
#[derive(Clone, Debug)]
pub struct HeavyLigand {
    /// The heavy-atom graph, atoms named `a0`, `a1`, ... in structure order.
    pub component: Component,
    /// Atom index in the structure of each component atom.
    pub atoms: Vec<usize>,
}

/// How the atoms of a model ligand were paired with the reference's.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LigandPairing {
    /// Every heavy atom name occurs once in both: paired by name.
    ByName,
    /// Names did not pair, but the heavy-atom elements agree in file order.
    ByOrder,
}

/// Why two ligands could not be compared.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
pub enum LigandCompareError {
    /// The reference ligand has no heavy atoms.
    #[error("the reference ligand has no heavy atoms")]
    NoHeavyAtoms,
    /// The two ligands have different numbers of heavy atoms.
    #[error("{reference} heavy atoms in the reference ligand, {model} in the prediction's")]
    AtomCount {
        /// Heavy atoms in the reference.
        reference: usize,
        /// Heavy atoms in the model.
        model: usize,
    },
    /// Neither names nor element order pair the atoms.
    #[error("the prediction's ligand atoms do not match the reference's")]
    Unpaired,
    /// The symmetry comparison failed (too many mappings, mismatched lengths).
    #[error("{0}")]
    Symmetry(String),
}

fn heavy_indices(ligand: &Structure) -> Vec<usize> {
    ligand
        .atoms()
        .into_iter()
        .filter(|a| a.element().is_some_and(|e| e.atomic_number() != 1))
        .map(|a| a.index().as_usize())
        .collect()
}

fn slots(heavy: &[usize]) -> HashMap<usize, usize> {
    heavy.iter().enumerate().map(|(s, a)| (*a, s)).collect()
}

/// The heavy atoms of `ligand` as a chemical component, or `None` when it has none.
#[must_use]
pub fn heavy_atom_component(ligand: &Structure) -> Option<HeavyLigand> {
    let heavy = heavy_indices(ligand);
    if heavy.is_empty() {
        return None;
    }
    let slot = slots(&heavy);
    let name = |s: usize| -> Box<str> { format!("a{s}").into() };
    let atoms: Vec<ComponentAtom> = heavy
        .iter()
        .enumerate()
        .filter_map(|(s, a)| {
            let element = ligand.atom_at(*a)?.element()?;
            Some(ComponentAtom {
                name: name(s),
                alternate_name: None,
                element,
                charge: 0,
                aromatic: false,
                leaving: false,
                stereo: None,
            })
        })
        .collect();
    let bonds: Vec<ComponentBond> = ligand
        .bonds()
        .iter()
        .filter_map(|b| {
            let (x, y) = (
                *slot.get(&b.atom_a.as_usize())?,
                *slot.get(&b.atom_b.as_usize())?,
            );
            Some(ComponentBond {
                atom_a: name(x),
                atom_b: name(y),
                order: b.order,
                aromatic: matches!(b.order, BondOrder::Aromatic),
                stereo: None,
            })
        })
        .collect();
    Some(HeavyLigand {
        component: Component {
            id: "ligand".into(),
            name: "ligand".into(),
            kind: ComponentKind::NonPolymer,
            parent: None,
            one_letter_code: None,
            formula: None,
            atoms: Arc::from(atoms),
            bonds: Arc::from(bonds),
            ideal_coordinates: None,
            model_coordinates: None,
        },
        atoms: heavy,
    })
}

/// The heavy-atom elements and bonds of a ligand as a comparable key: two poses of one ligand
/// have the same key.
#[must_use]
pub fn heavy_graph_key(ligand: &Structure) -> (Vec<u8>, Vec<(usize, usize)>) {
    let heavy = heavy_indices(ligand);
    let slot = slots(&heavy);
    let elements = heavy
        .iter()
        .filter_map(|i| {
            ligand
                .atom_at(*i)
                .and_then(molframe_core::structure::AtomRef::element)
                .map(molframe_core::Element::atomic_number)
        })
        .collect();
    let mut bonds: Vec<(usize, usize)> = ligand
        .bonds()
        .iter()
        .filter_map(|b| {
            let (x, y) = (
                *slot.get(&b.atom_a.as_usize())?,
                *slot.get(&b.atom_b.as_usize())?,
            );
            Some((x.min(y), x.max(y)))
        })
        .collect();
    bonds.sort_unstable();
    (elements, bonds)
}

fn coordinates_of(ligand: &Structure, heavy: &[usize]) -> Vec<[f32; 3]> {
    let all = ligand.coordinates();
    heavy.iter().filter_map(|i| all.get(*i).copied()).collect()
}

/// RMSD between two poses of one ligand, heavy atoms only, with no superposition and the best of
/// the ligand's symmetry-equivalent atom mappings (at most `mapping_limit` tried). `None` when
/// they are not the same ligand (a different heavy-atom graph) or the comparison failed.
#[must_use]
pub fn pose_rmsd(a: &Structure, b: &Structure, mapping_limit: usize) -> Option<f64> {
    if heavy_graph_key(a) != heavy_graph_key(b) {
        return None;
    }
    let first = heavy_atom_component(a)?;
    let heavy_b = heavy_indices(b);
    if first.atoms.len() != heavy_b.len() {
        return None;
    }
    let (ra, rb) = (coordinates_of(a, &first.atoms), coordinates_of(b, &heavy_b));
    ligand_symmetry_rmsd(&ra, &rb, &first.component, mapping_limit)
        .ok()
        .map(|r| r.rmsd)
}

/// Heavy atoms of `ligand` as (name, atomic number, position).
fn named_heavy_atoms(ligand: &Structure) -> Vec<(String, u8, [f32; 3])> {
    ligand
        .atoms()
        .into_iter()
        .filter_map(|a| {
            let e = a.element().filter(|e| e.atomic_number() > 1)?;
            Some((
                a.name()?.trim().to_owned(),
                e.atomic_number(),
                a.position()?,
            ))
        })
        .collect()
}

/// Applies a column-major 4x4 placement to a point.
fn place(m: &[f32; 16], p: [f32; 3]) -> [f32; 3] {
    std::array::from_fn(|i| m[i] * p[0] + m[4 + i] * p[1] + m[8 + i] * p[2] + m[12 + i])
}

/// Heavy-atom RMSD of `model` against `reference` after `model` is carried by `placement`
/// (column-major 4x4), with the best of the reference ligand's symmetry-equivalent mappings. The
/// model needs no bonds and no matching atom order: its atoms are paired with the reference's by
/// name, else by order when the elements agree.
///
/// # Errors
/// Why the two could not be compared (different atoms, too many symmetries).
pub fn placed_rmsd(
    reference: &Structure,
    model: &Structure,
    placement: &[f32; 16],
    mapping_limit: usize,
) -> Result<(f64, LigandPairing), LigandCompareError> {
    let heavy = heavy_atom_component(reference).ok_or(LigandCompareError::NoHeavyAtoms)?;
    let refs = coordinates_of(reference, &heavy.atoms);
    let ref_atoms = named_heavy_atoms(reference);
    let mod_atoms = named_heavy_atoms(model);
    if ref_atoms.len() != mod_atoms.len() || ref_atoms.len() != refs.len() {
        return Err(LigandCompareError::AtomCount {
            reference: ref_atoms.len(),
            model: mod_atoms.len(),
        });
    }
    let unique = |v: &[(String, u8, [f32; 3])]| {
        let mut n: Vec<&str> = v.iter().map(|a| a.0.as_str()).collect();
        n.sort_unstable();
        n.windows(2).all(|w| w[0] != w[1])
    };
    let by_name: Option<Vec<[f32; 3]>> = (unique(&ref_atoms) && unique(&mod_atoms))
        .then(|| {
            ref_atoms
                .iter()
                .map(|r| {
                    mod_atoms
                        .iter()
                        .find(|m| m.0 == r.0 && m.1 == r.1)
                        .map(|m| m.2)
                })
                .collect()
        })
        .flatten();
    let (coords, pairing) = match by_name {
        Some(c) => (c, LigandPairing::ByName),
        None if ref_atoms.iter().zip(&mod_atoms).all(|(r, m)| r.1 == m.1) => (
            mod_atoms.iter().map(|m| m.2).collect(),
            LigandPairing::ByOrder,
        ),
        None => return Err(LigandCompareError::Unpaired),
    };
    let moved: Vec<[f32; 3]> = coords.iter().map(|p| place(placement, *p)).collect();
    ligand_symmetry_rmsd(&refs, &moved, &heavy.component, mapping_limit)
        .map(|r| (r.rmsd, pairing))
        .map_err(|e| LigandCompareError::Symmetry(e.to_string()))
}

#[cfg(test)]
#[path = "ligand_tests.rs"]
mod tests;
