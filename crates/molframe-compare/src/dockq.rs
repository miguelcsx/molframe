//! `DockQ`, the published quality score for a two-body docking model.
//!
//! This follows Basu and Wallner (2016). Three views of how well a model
//! reproduces a native complex are combined:
//!
//! - `fnat`, the fraction of native residue-residue contacts the model keeps. Two
//!   residues of different chains are in contact when any heavy atom of one lies
//!   closer than the contact distance (5 Å) to any heavy atom of the other. Each
//!   residue pair counts once however many atom pairs touch.
//! - `LRMSD`, the backbone (N, CA, C, O) RMSD of the ligand after the model's
//!   receptor backbone is superposed on the native receptor.
//! - `iRMSD`, the backbone RMSD over the interface residues after superposing on
//!   those. Interface residues are the native residues having any heavy atom
//!   closer than the interface distance (10 Å) to a heavy atom of the other chain.
//!
//! ```text
//! DockQ = (fnat + 1 / (1 + (LRMSD / 8.5)^2) + 1 / (1 + (iRMSD / 1.5)^2)) / 3
//! ```
//!
//! The ligand is the chain the caller names as such; the published program picks
//! the smaller chain, and callers wanting that should name it. A native without a
//! single residue contact has no `fnat` and is refused rather than scored.
//!
//! Model and native must share atom numbering. The atom-contact variant this
//! crate used to call `DockQ` lives in [`atom_contact_docking_score`] and is not
//! comparable with published values.

use std::collections::BTreeSet;

use molframe_core::contract::Namespace;
use molframe_core::structure::Structure;
use molframe_geom::{rmsd, superpose};

use crate::failure::CompareError;
use crate::interface::{AtomInfo, UNANNOTATED, atom_infos, chain_atoms};
use crate::numeric::usize_to_f64;

#[path = "dockq_atom_contact.rs"]
mod atom_contact;

pub use atom_contact::{
    AtomContactDockingOptions, AtomContactDockingScore, atom_contact_docking_score,
    atom_contact_docking_score_in_namespace,
};

/// Explicit definition of a `DockQ` calculation.
///
/// [`DockQOptions::published`] gives the published definition. The two distance
/// cutoffs are independent: a residue pair can be an interface neighbour without
/// being a contact.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DockQOptions {
    /// Heavy-atom distance below which two residues are in contact, in angstroms.
    pub contact_distance: f32,
    /// Heavy-atom distance below which a residue belongs to the interface, in
    /// angstroms.
    pub interface_distance: f32,
    /// Length scale used to squash ligand RMSD, in angstroms.
    pub ligand_scale: f64,
    /// Length scale used to squash interface RMSD, in angstroms.
    pub interface_scale: f64,
}

impl DockQOptions {
    /// The published parameters: 5 Å contacts, a 10 Å interface and length
    /// scales of 8.5 Å and 1.5 Å.
    #[must_use]
    pub const fn published() -> Self {
        Self {
            contact_distance: 5.0,
            interface_distance: 10.0,
            ligand_scale: 8.5,
            interface_scale: 1.5,
        }
    }

    fn validate(self) -> Result<Self, CompareError> {
        let cutoffs = [self.contact_distance, self.interface_distance];
        if cutoffs
            .iter()
            .any(|cutoff| !cutoff.is_finite() || *cutoff <= 0.0)
        {
            return Err(CompareError::InvalidDistanceCutoff);
        }
        if !self.ligand_scale.is_finite()
            || self.ligand_scale <= 0.0
            || !self.interface_scale.is_finite()
            || self.interface_scale <= 0.0
        {
            return Err(CompareError::InvalidLengthScale);
        }
        Ok(self)
    }
}

/// A `DockQ` score and the three components it averages.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DockQ {
    /// Fraction of native residue-residue contacts kept by the model.
    pub fnat: f64,
    /// Backbone RMSD of the ligand after superposing the receptor backbone.
    pub ligand_rmsd: f64,
    /// Backbone RMSD over the interface residues.
    pub interface_rmsd: f64,
    /// The combined score in `[0, 1]`.
    pub score: f64,
}

/// Scores a docking `model` against the `native` complex.
///
/// `receptor` and `ligand` name the two chains; the receptor is the frame the
/// ligand is measured against. Identical structures score `1.0`.
///
/// # Errors
///
/// Returns [`CompareError::LengthMismatch`] when the structures disagree on atom
/// count, [`CompareError::InvalidDistanceCutoff`] or
/// [`CompareError::InvalidLengthScale`] when `options` are not meaningful,
/// [`CompareError::NoComparablePairs`] when the native has no residue contact or
/// a chain has no backbone atoms, and [`CompareError::Superpose`] when a
/// superposition has too few atoms.
pub fn dockq(
    model: &Structure,
    native: &Structure,
    receptor: &str,
    ligand: &str,
    options: DockQOptions,
) -> Result<DockQ, CompareError> {
    dockq_in_namespace(model, native, receptor, ligand, Namespace::Label, options)
}

/// Scores a docking model with chain names interpreted in one namespace.
///
/// # Errors
///
/// Returns the same errors as [`dockq`] and rejects an explicit namespace.
pub fn dockq_in_namespace(
    model: &Structure,
    native: &Structure,
    receptor: &str,
    ligand: &str,
    namespace: Namespace,
    options: DockQOptions,
) -> Result<DockQ, CompareError> {
    if model.atom_count() != native.atom_count() {
        return Err(CompareError::LengthMismatch {
            model: model.atom_count() as usize,
            reference: native.atom_count() as usize,
        });
    }
    let receptor_atoms = chain_atoms(native, receptor, namespace)?;
    let ligand_atoms = chain_atoms(native, ligand, namespace)?;
    dockq_on_positions(
        model.positions(),
        native.positions(),
        &receptor_atoms,
        &ligand_atoms,
        &atom_infos(native),
        options,
    )
}

/// Scores model coordinates against native ones, row for row.
///
/// `receptor_rows` and `ligand_rows` index rows of both coordinate sets and
/// `info` annotates each row; a row beyond `info` takes part in nothing.
///
/// # Errors
///
/// Returns the errors of [`dockq`] for invalid options or too few atoms.
pub(crate) fn dockq_on_positions(
    model: &[[f32; 3]],
    native: &[[f32; 3]],
    receptor_rows: &[usize],
    ligand_rows: &[usize],
    info: &[AtomInfo],
    options: DockQOptions,
) -> Result<DockQ, CompareError> {
    let options = options.validate()?;
    let sides = Sides {
        receptor: receptor_rows,
        ligand: ligand_rows,
        info,
    };
    let native_contacts = sides.residue_pairs(native, options.contact_distance);
    if native_contacts.is_empty() {
        return Err(CompareError::NoComparablePairs);
    }
    let model_contacts = sides.residue_pairs(model, options.contact_distance);
    let kept = native_contacts.intersection(&model_contacts).count();
    let fnat = usize_to_f64(kept) / usize_to_f64(native_contacts.len());

    let ligand_rmsd = sides.ligand_rmsd(model, native)?;
    let interface = sides.interface_residues(native, options.interface_distance);
    let interface_rmsd = sides.interface_rmsd(model, native, &interface)?;

    let score = (fnat
        + squash(ligand_rmsd, options.ligand_scale)
        + squash(interface_rmsd, options.interface_scale))
        / 3.0;
    Ok(DockQ {
        fnat,
        ligand_rmsd,
        interface_rmsd,
        score,
    })
}

/// The two partners and the annotations of their rows.
struct Sides<'a> {
    receptor: &'a [usize],
    ligand: &'a [usize],
    info: &'a [AtomInfo],
}

impl Sides<'_> {
    fn info(&self, row: usize) -> AtomInfo {
        match self.info.get(row) {
            Some(&found) => found,
            None => UNANNOTATED,
        }
    }

    /// Heavy atoms of one partner as `(residue, position)`.
    fn heavy(&self, rows: &[usize], positions: &[[f32; 3]]) -> Vec<(usize, [f32; 3])> {
        rows.iter()
            .filter_map(|&row| {
                let info = self.info(row);
                let point = positions.get(row)?;
                info.heavy.then_some((info.residue, *point))
            })
            .collect()
    }

    /// Residue pairs `(receptor residue, ligand residue)` having any heavy-atom
    /// pair strictly closer than `cutoff`.
    fn residue_pairs(&self, positions: &[[f32; 3]], cutoff: f32) -> BTreeSet<(usize, usize)> {
        let limit = f64::from(cutoff) * f64::from(cutoff);
        let ligand = self.heavy(self.ligand, positions);
        let mut pairs = BTreeSet::new();
        for (r, a) in self.heavy(self.receptor, positions) {
            for &(l, b) in &ligand {
                if molframe_geom::distance_squared(a, b) < limit {
                    pairs.insert((r, l));
                }
            }
        }
        pairs
    }

    /// Residues of either partner within `cutoff` of the other partner.
    fn interface_residues(&self, native: &[[f32; 3]], cutoff: f32) -> BTreeSet<usize> {
        let mut residues = BTreeSet::new();
        for (r, l) in self.residue_pairs(native, cutoff) {
            residues.insert(r);
            residues.insert(l);
        }
        residues
    }

    /// Backbone rows among `rows`, optionally limited to some residues.
    fn backbone(&self, rows: &[usize], only: Option<&BTreeSet<usize>>) -> Vec<usize> {
        rows.iter()
            .copied()
            .filter(|&row| {
                let info = self.info(row);
                info.backbone && only.is_none_or(|set| set.contains(&info.residue))
            })
            .collect()
    }

    /// Ligand backbone RMSD after superposing the receptor backbone.
    fn ligand_rmsd(&self, model: &[[f32; 3]], native: &[[f32; 3]]) -> Result<f64, CompareError> {
        let receptor = self.backbone(self.receptor, None);
        let ligand = self.backbone(self.ligand, None);
        if receptor.is_empty() || ligand.is_empty() {
            return Err(CompareError::NoComparablePairs);
        }
        let fit = superpose(&gather(model, &receptor), &gather(native, &receptor))
            .map_err(CompareError::Superpose)?;
        let moved: Vec<[f32; 3]> = gather(model, &ligand)
            .into_iter()
            .map(|point| fit.transform.apply(point))
            .collect();
        rmsd(&moved, &gather(native, &ligand)).map_err(CompareError::Superpose)
    }

    /// Backbone RMSD of the interface residues after superposing on them.
    fn interface_rmsd(
        &self,
        model: &[[f32; 3]],
        native: &[[f32; 3]],
        interface: &BTreeSet<usize>,
    ) -> Result<f64, CompareError> {
        let mut rows = self.backbone(self.receptor, Some(interface));
        rows.extend(self.backbone(self.ligand, Some(interface)));
        if rows.is_empty() {
            return Err(CompareError::NoComparablePairs);
        }
        let fit = superpose(&gather(model, &rows), &gather(native, &rows))
            .map_err(CompareError::Superpose)?;
        Ok(fit.rmsd)
    }
}

/// Squashes a distance onto `(0, 1]` with a length scale.
fn squash(distance: f64, scale: f64) -> f64 {
    let ratio = distance / scale;
    1.0 / (1.0 + ratio * ratio)
}

/// Collects the positions at the given indices, in order.
fn gather(positions: &[[f32; 3]], indices: &[usize]) -> Vec<[f32; 3]> {
    indices
        .iter()
        .filter_map(|&index| positions.get(index).copied())
        .collect()
}

#[cfg(test)]
#[path = "dockq_tests.rs"]
mod tests;
