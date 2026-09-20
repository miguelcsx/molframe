//! Structure-first analysis operations, as methods.
//!
//! Every method here forwards to a function in [`crate::analysis`]; read that
//! module for what each operation computes. What the method form adds is the
//! two defaults a convenience call wants and the kernel cannot assume: the
//! backend is [`SpatialBackend::Auto`], which picks from the workload's shape,
//! and the context is [`ExecutionContext::default`], the bounded profile. The
//! `_with` form takes both explicitly and is the one a caller who cares about
//! either should write.

use crate::ExecutionContext;
use crate::structure::{Selection, Structure};
use molframe_chem::ComponentProvider;
use molframe_spatial::{SpatialBackend, SpatialError};

/// Analysis operations over a structure.
///
/// # Examples
///
/// ```
/// # #[cfg(feature = "analysis")]
/// # {
/// use molframe::prelude::*;
///
/// # const PDB: &str = "\
/// # ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N
/// # ATOM      2  CA  ALA A   1      12.560   6.195  -6.504  1.00  0.00           C
/// # ATOM      3  C   ALA A   1      13.100   7.520  -6.504  1.00  0.00           C
/// # END
/// # ";
/// let (structure, _) = read_bytes(PDB.into(), Some("1abc.pdb"), &ReadOptions::new())?;
/// let contacts = structure.atom_contacts(4.5)?;
/// assert!(!contacts.is_empty());
/// # }
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub trait AnalysisExt {
    /// Every unordered atom pair within `cutoff` ångström.
    ///
    /// # Errors
    ///
    /// Returns a spatial error for an invalid cutoff or a backend that cannot
    /// serve the structure.
    fn atom_contacts(&self, cutoff: f32) -> Result<molframe_analysis::ContactTable, SpatialError>;

    /// [`Self::atom_contacts`] with the backend and context given explicitly.
    ///
    /// # Errors
    ///
    /// As [`Self::atom_contacts`].
    fn atom_contacts_with(
        &self,
        cutoff: f32,
        backend: SpatialBackend,
        context: &ExecutionContext,
    ) -> Result<molframe_analysis::ContactTable, SpatialError>;

    /// Every pair across two selections within `cutoff` ångström.
    ///
    /// # Errors
    ///
    /// As [`Self::atom_contacts`].
    fn atom_contacts_between(
        &self,
        left: &Selection,
        right: &Selection,
        cutoff: f32,
    ) -> Result<molframe_analysis::ContactTable, SpatialError>;

    /// [`Self::atom_contacts_between`] with the backend and context given
    /// explicitly.
    ///
    /// # Errors
    ///
    /// As [`Self::atom_contacts`].
    fn atom_contacts_between_with(
        &self,
        left: &Selection,
        right: &Selection,
        cutoff: f32,
        backend: SpatialBackend,
        context: &ExecutionContext,
    ) -> Result<molframe_analysis::ContactTable, SpatialError>;

    /// Streams every unordered atom pair within `cutoff` to `emit`.
    ///
    /// # Errors
    ///
    /// As [`Self::atom_contacts`].
    fn visit_atom_contacts<F: FnMut(molframe_analysis::Contact)>(
        &self,
        cutoff: f32,
        emit: F,
    ) -> Result<(), SpatialError>;

    /// [`Self::visit_atom_contacts`] with the backend and context given
    /// explicitly.
    ///
    /// # Errors
    ///
    /// As [`Self::atom_contacts`].
    fn visit_atom_contacts_with<F: FnMut(molframe_analysis::Contact)>(
        &self,
        cutoff: f32,
        backend: SpatialBackend,
        context: &ExecutionContext,
        emit: F,
    ) -> Result<(), SpatialError>;

    /// The residue-level contact map within `cutoff` ångström.
    ///
    /// # Errors
    ///
    /// As [`Self::atom_contacts`].
    fn residue_contact_map(
        &self,
        cutoff: f32,
        min_separation: u32,
    ) -> Result<molframe_analysis::ContactMap, SpatialError>;

    /// [`Self::residue_contact_map`] with the backend and context given
    /// explicitly.
    ///
    /// # Errors
    ///
    /// As [`Self::atom_contacts`].
    fn residue_contact_map_with(
        &self,
        cutoff: f32,
        min_separation: u32,
        backend: SpatialBackend,
        context: &ExecutionContext,
    ) -> Result<molframe_analysis::ContactMap, SpatialError>;

    /// Hydrogen bonds under `options`.
    ///
    /// # Errors
    ///
    /// Returns a hydrogen-bond error for an invalid geometry threshold or a
    /// structure the chemistry cannot classify.
    fn hydrogen_bonds(
        &self,
        options: molframe_analysis::HydrogenBondOptions,
    ) -> Result<molframe_analysis::HydrogenBondTable, molframe_analysis::HydrogenBondError>;

    /// [`Self::hydrogen_bonds`] with the context given explicitly.
    ///
    /// # Errors
    ///
    /// As [`Self::hydrogen_bonds`].
    fn hydrogen_bonds_with(
        &self,
        options: molframe_analysis::HydrogenBondOptions,
        context: &ExecutionContext,
    ) -> Result<molframe_analysis::HydrogenBondTable, molframe_analysis::HydrogenBondError>;

    /// Water bridges under `options`.
    ///
    /// # Errors
    ///
    /// As [`Self::hydrogen_bonds`].
    fn water_bridges(
        &self,
        options: molframe_analysis::WaterBridgeOptions,
    ) -> Result<molframe_analysis::WaterBridgeTable, molframe_analysis::HydrogenBondError>;

    /// [`Self::water_bridges`] with the context given explicitly.
    ///
    /// # Errors
    ///
    /// As [`Self::hydrogen_bonds`].
    fn water_bridges_with(
        &self,
        options: molframe_analysis::WaterBridgeOptions,
        context: &ExecutionContext,
    ) -> Result<molframe_analysis::WaterBridgeTable, molframe_analysis::HydrogenBondError>;

    /// Salt bridges within `max_distance` ångström.
    ///
    /// # Errors
    ///
    /// As [`Self::atom_contacts`].
    fn salt_bridges(
        &self,
        max_distance: f32,
    ) -> Result<molframe_analysis::SaltBridgeTable, SpatialError>;

    /// [`Self::salt_bridges`] with the backend and context given explicitly.
    ///
    /// # Errors
    ///
    /// As [`Self::atom_contacts`].
    fn salt_bridges_with(
        &self,
        max_distance: f32,
        backend: SpatialBackend,
        context: &ExecutionContext,
    ) -> Result<molframe_analysis::SaltBridgeTable, SpatialError>;

    /// π-stacking interactions under `options`.
    ///
    /// # Errors
    ///
    /// Returns a stacking error for invalid ring-geometry thresholds.
    fn pi_stacking(
        &self,
        options: molframe_analysis::PiStackingOptions,
    ) -> Result<molframe_analysis::PiStackingTable, molframe_analysis::PiStackingError>;

    /// Cation–π interactions under `options`.
    ///
    /// # Errors
    ///
    /// Returns a cation–π error for invalid geometry thresholds.
    fn cation_pi(
        &self,
        options: molframe_analysis::CationPiOptions,
    ) -> Result<molframe_analysis::CationPiTable, molframe_analysis::CationPiError>;

    /// Base pairs under `options`, with `provider` supplying the chemistry.
    ///
    /// # Errors
    ///
    /// Returns a base-pair error for a component the provider does not define.
    fn base_pairs(
        &self,
        provider: &dyn ComponentProvider,
        options: molframe_analysis::BasePairOptions,
    ) -> Result<Vec<molframe_analysis::BasePair>, molframe_analysis::BasePairError>;

    /// [`Self::base_pairs`] with the context given explicitly.
    ///
    /// # Errors
    ///
    /// As [`Self::base_pairs`].
    fn base_pairs_with(
        &self,
        provider: &dyn ComponentProvider,
        options: molframe_analysis::BasePairOptions,
        context: &ExecutionContext,
    ) -> Result<Vec<molframe_analysis::BasePair>, molframe_analysis::BasePairError>;

    /// Secondary structure assignment under `options`.
    ///
    /// # Errors
    ///
    /// Returns an assignment error when the structure has no polymer trace the
    /// algorithm can follow.
    fn secondary_structure(
        &self,
        options: &molframe_analysis::DsspOptions,
    ) -> Result<molframe_analysis::SseTable, molframe_analysis::DsspError>;

    /// The residues at the interface of two chains, within `cutoff` ångström.
    ///
    /// # Errors
    ///
    /// As [`Self::atom_contacts`].
    fn chain_interface(
        &self,
        first: &str,
        second: &str,
        cutoff: f32,
    ) -> Result<Vec<molframe_core::index::ResidueIndex>, SpatialError>;

    /// [`Self::chain_interface`] with the backend and context given
    /// explicitly.
    ///
    /// # Errors
    ///
    /// As [`Self::atom_contacts`].
    fn chain_interface_with(
        &self,
        first: &str,
        second: &str,
        cutoff: f32,
        backend: SpatialBackend,
        context: &ExecutionContext,
    ) -> Result<Vec<molframe_core::index::ResidueIndex>, SpatialError>;

    /// The fraction of `target`'s native contacts this structure preserves.
    ///
    /// The receiver is the reference, as it is for the free function.
    ///
    /// # Errors
    ///
    /// As [`Self::atom_contacts`].
    fn native_contact_fraction(
        &self,
        target: &Structure,
        cutoff: f32,
        tolerance: f32,
    ) -> Result<molframe_analysis::NativeContacts, molframe_analysis::NativeError>;

    /// [`Self::native_contact_fraction`] with the backend and context given
    /// explicitly.
    ///
    /// # Errors
    ///
    /// As [`Self::atom_contacts`].
    fn native_contact_fraction_with(
        &self,
        target: &Structure,
        cutoff: f32,
        tolerance: f32,
        backend: SpatialBackend,
        context: &ExecutionContext,
    ) -> Result<molframe_analysis::NativeContacts, molframe_analysis::NativeError>;

    /// The half-sphere exposure of every atom, at `radius` ångström.
    ///
    /// # Errors
    ///
    /// Returns an exposure error for an invalid radius.
    fn half_sphere_exposure(
        &self,
        radius: f32,
    ) -> Result<Vec<molframe_analysis::HalfSphereExposure>, molframe_analysis::HseError>;

    /// [`Self::half_sphere_exposure`] with the backend and context given
    /// explicitly.
    ///
    /// # Errors
    ///
    /// As [`Self::half_sphere_exposure`].
    fn half_sphere_exposure_with(
        &self,
        radius: f32,
        backend: SpatialBackend,
        context: &ExecutionContext,
    ) -> Result<Vec<molframe_analysis::HalfSphereExposure>, molframe_analysis::HseError>;

    /// The backbone torsions of every nucleotide.
    ///
    /// # Errors
    ///
    /// Returns a torsion error when a nucleotide's trace is incomplete.
    fn nucleic_torsions(
        &self,
    ) -> Result<Vec<molframe_analysis::NucleicTorsions>, molframe_analysis::NucleicTorsionError>;

    /// Solvent-accessible contacts over caller-supplied per-atom `radii`.
    ///
    /// # Errors
    ///
    /// Returns a surface error for an invalid probe, radius or sampling
    /// density.
    fn surface_contacts(
        &self,
        radii: &[f32],
        options: molframe_analysis::SurfaceContactOptions,
    ) -> Result<Vec<molframe_analysis::Contact>, molframe_surface::SasaError>;

    /// [`Self::surface_contacts`] with the context given explicitly.
    ///
    /// # Errors
    ///
    /// As [`Self::surface_contacts`].
    fn surface_contacts_with(
        &self,
        radii: &[f32],
        options: molframe_analysis::SurfaceContactOptions,
        context: &ExecutionContext,
    ) -> Result<Vec<molframe_analysis::Contact>, molframe_surface::SasaError>;
}
