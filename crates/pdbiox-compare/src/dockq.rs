//! `DockQ`, a single quality score for a two-body docking model.
//!
//! `DockQ` blends three views of how well a model reproduces a native complex: the
//! fraction of native interface contacts it keeps (Fnat), how far its ligand
//! sits once the receptor is superposed (the ligand RMSD), and how well the
//! interface itself superposes (the interface RMSD). The two distances are
//! squashed onto `(0, 1]` by a length scale and averaged with Fnat, giving a
//! score that is 1 for a perfect model and near 0 for an unrelated one.
//!
//! Model and native must share atom numbering — the model is a pose of the same
//! atoms — and the two chains are named. Cost is one receptor superposition, one
//! interface superposition, and a pass over the interface contacts.

use pdbiox_core::contract::Namespace;
use pdbiox_core::structure::Structure;
use pdbiox_geom::{rmsd, superpose};

use crate::failure::CompareError;
use crate::interface::{chain_atoms, contacts};
use crate::numeric::{u64_to_f64, usize_to_f64};

/// Explicit definition of a `DockQ` calculation.
///
/// The library deliberately supplies no implicit literature or compatibility
/// defaults. Callers choose the contact definition and both score scales, and
/// can therefore record them alongside the result.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DockQOptions {
    /// Maximum native inter-chain distance counted as a contact, in angstroms.
    pub contact_distance: f32,
    /// Length scale used to squash ligand RMSD, in angstroms.
    pub ligand_scale: f64,
    /// Length scale used to squash interface RMSD, in angstroms.
    pub interface_scale: f64,
}

impl DockQOptions {
    fn validate(self) -> Result<Self, CompareError> {
        if !self.contact_distance.is_finite() || self.contact_distance <= 0.0 {
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
    /// Fraction of native interface contacts kept by the model.
    pub fnat: f64,
    /// Ligand RMSD after superposing the receptor.
    pub ligand_rmsd: f64,
    /// Interface RMSD.
    pub interface_rmsd: f64,
    /// The combined score in `(0, 1]`.
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
/// [`CompareError::InvalidLengthScale`] when `options` are not physically
/// meaningful, and [`CompareError::Superpose`] when a chain or interface has too
/// few atoms to superpose.
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
    let options = options.validate()?;
    if model.atom_count() != native.atom_count() {
        return Err(CompareError::LengthMismatch {
            model: model.atom_count() as usize,
            reference: native.atom_count() as usize,
        });
    }
    let receptor_atoms = chain_atoms(native, receptor, namespace)?;
    let ligand_atoms = chain_atoms(native, ligand, namespace)?;
    let model_positions = model.positions();
    let native_positions = native.positions();

    let native_contacts = contacts(
        native_positions,
        &receptor_atoms,
        &ligand_atoms,
        options.contact_distance,
    );
    let fnat = fraction_kept(model_positions, &native_contacts, options.contact_distance);

    let ligand_rmsd = ligand_rmsd(
        model_positions,
        native_positions,
        &receptor_atoms,
        &ligand_atoms,
    )?;
    let interface_rmsd = interface_rmsd(model_positions, native_positions, &native_contacts)?;

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

/// Squashes a distance onto `(0, 1]` with a length scale.
fn squash(distance: f64, scale: f64) -> f64 {
    let ratio = distance / scale;
    1.0 / (1.0 + ratio * ratio)
}

/// The fraction of native contacts still within the cutoff in the model.
fn fraction_kept(
    positions: &[[f32; 3]],
    native_contacts: &[(usize, usize)],
    contact_distance: f32,
) -> f64 {
    if native_contacts.is_empty() {
        return 1.0;
    }
    let cutoff_squared = f64::from(contact_distance) * f64::from(contact_distance);
    let mut kept = 0u64;
    for &(r, l) in native_contacts {
        let (Some(&a), Some(&b)) = (positions.get(r), positions.get(l)) else {
            continue;
        };
        if pdbiox_geom::distance_squared(a, b) <= cutoff_squared {
            kept += 1;
        }
    }
    u64_to_f64(kept) / usize_to_f64(native_contacts.len())
}

/// The ligand RMSD after superposing the model receptor onto the native.
fn ligand_rmsd(
    model: &[[f32; 3]],
    native: &[[f32; 3]],
    receptor: &[usize],
    ligand: &[usize],
) -> Result<f64, CompareError> {
    let model_receptor = gather(model, receptor);
    let native_receptor = gather(native, receptor);
    let fit = superpose(&model_receptor, &native_receptor).map_err(CompareError::Superpose)?;

    let moved: Vec<[f32; 3]> = ligand
        .iter()
        .filter_map(|&index| model.get(index).map(|&point| fit.transform.apply(point)))
        .collect();
    let native_ligand = gather(native, ligand);
    rmsd(&moved, &native_ligand).map_err(CompareError::Superpose)
}

/// The interface RMSD: superposition of the atoms in the native contacts.
fn interface_rmsd(
    model: &[[f32; 3]],
    native: &[[f32; 3]],
    native_contacts: &[(usize, usize)],
) -> Result<f64, CompareError> {
    let mut interface: Vec<usize> = Vec::new();
    for &(r, l) in native_contacts {
        interface.push(r);
        interface.push(l);
    }
    interface.sort_unstable();
    interface.dedup();

    let model_interface = gather(model, &interface);
    let native_interface = gather(native, &interface);
    let fit = superpose(&model_interface, &native_interface).map_err(CompareError::Superpose)?;
    Ok(fit.rmsd)
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
