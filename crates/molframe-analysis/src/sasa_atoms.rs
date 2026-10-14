//! Which atoms of a structure take part in a solvent-accessible calculation.
//!
//! One rule shared by the per-residue areas and the streamed totals: the atoms
//! the heteroatom options include, hydrogens left out, one alternate conformer
//! per residue (the first label met), and a radius from the options' table.

use molframe_chem::included_atoms;
use molframe_core::structure::Structure;
use molframe_surface::SasaOptions;

/// The radius of every atom, `NaN` for an atom that is left out.
pub(crate) struct MeasuredRadii {
    /// One radius per atom of the structure.
    pub(crate) radii: Vec<f32>,
    /// Included atoms (not hydrogens, not other conformers) with no radius.
    pub(crate) missing: usize,
}

/// The radii of the atoms `options` measure, in atom order.
pub(crate) fn measured_radii(structure: &Structure, options: &SasaOptions) -> MeasuredRadii {
    let mut radii = options.atom_radii(structure);
    let included = included_atoms(structure, options.hetero);
    let mut missing = 0;
    for chain in structure.data().chains() {
        for residue in chain.residues() {
            let mut conformer: Option<&str> = None;
            for atom in residue.atoms() {
                let index = atom.index().as_usize();
                if !matches!(included.get(index), Some(true)) {
                    continue;
                }
                let hydrogen = atom
                    .element()
                    .is_some_and(|element| element.atomic_number() <= 1);
                let other_conformer = atom
                    .alt_label()
                    .is_some_and(|alt| *conformer.get_or_insert(alt) != alt);
                let Some(slot) = radii.get_mut(index) else {
                    continue;
                };
                if hydrogen || other_conformer {
                    *slot = f32::NAN;
                } else if !(slot.is_finite() && *slot > 0.0) {
                    *slot = f32::NAN;
                    missing += 1;
                }
            }
        }
    }
    MeasuredRadii { radii, missing }
}
