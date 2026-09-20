//! Structure-first validation operations, as methods.
//!
//! Every method here forwards to a function in [`crate::validation`]. Only
//! [`ValidationExt::clashes`] has an `_with` twin, because it is the only
//! kernel here that both schedules work through a spatial backend and takes a
//! context. The others take the whole of their configuration as an options
//! value already, so a convenience form would have nothing to default.

use crate::structure::{Selection, Structure};
use crate::{AnalysisPolicy, ExecutionContext, Namespace};
use molframe_chem::{ComponentProvider, PolymerRoleProfile, RadiusSet};
use molframe_spatial::SpatialBackend;

/// Validation operations over a structure.
///
/// # Examples
///
/// ```
/// # #[cfg(feature = "validation")]
/// # {
/// use molframe::prelude::*;
///
/// # const PDB: &str = "\
/// # ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N
/// # ATOM      2  CA  ALA A   1      12.560   6.195  -6.504  1.00  0.00           C
/// # END
/// # ";
/// let (structure, _) = read_bytes(PDB.into(), Some("1abc.pdb"), &ReadOptions::new())?;
/// let flags = structure.quality_flags();
/// assert!(flags.len() <= structure.atom_count() as usize);
/// # }
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub trait ValidationExt {
    /// Steric clashes under `tolerance` ångström, against `radius_set`.
    ///
    /// # Errors
    ///
    /// Returns a spatial error for an invalid tolerance or a backend that
    /// cannot serve the structure.
    fn clashes(
        &self,
        tolerance: f32,
        radius_set: RadiusSet,
    ) -> Result<molframe_validate::ClashTable, molframe_spatial::SpatialError>;

    /// [`Self::clashes`] with the backend and context given explicitly.
    ///
    /// # Errors
    ///
    /// As [`Self::clashes`].
    fn clashes_with(
        &self,
        tolerance: f32,
        radius_set: RadiusSet,
        backend: SpatialBackend,
        context: &ExecutionContext,
    ) -> Result<molframe_validate::ClashTable, molframe_spatial::SpatialError>;

    /// Bond-length deviations beyond `tolerance` ångström.
    #[must_use]
    fn bond_length_deviations(&self, tolerance: f32) -> Vec<molframe_validate::BondDeviation>;

    /// The ligand geometry report at `tolerance` ångström.
    #[must_use]
    fn ligand_geometry(&self, tolerance: f32) -> molframe_validate::LigandGeometryReport;

    /// The ligand bond deviations beyond `tolerance` ångström.
    #[must_use]
    fn ligand_geometry_outliers(&self, tolerance: f32) -> Vec<molframe_validate::BondDeviation>;

    /// One quality flag per atom: occupancy and B-factor sanity.
    #[must_use]
    fn quality_flags(&self) -> Vec<molframe_validate::QualityFlag>;

    /// Atoms whose bond order sum exceeds their element's valence.
    #[must_use]
    fn overvalent_atoms(&self) -> Vec<molframe_validate::ValenceError>;

    /// Cis-peptide bonds whose ω deviates past `threshold_degrees`.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic when the structure has no polymer trace to walk.
    fn cis_peptides(
        &self,
        threshold_degrees: f64,
    ) -> Result<Vec<molframe_validate::CisPeptide>, molframe_core::diagnostic::Diagnostic>;

    /// Aromatic rings that are not planar under `options`.
    ///
    /// # Errors
    ///
    /// Returns a planarity error for invalid options.
    fn nonplanar_aromatic_rings(
        &self,
        options: molframe_validate::PlanarityOptions,
    ) -> Result<Vec<molframe_validate::PlanarityFlag>, molframe_validate::PlanarityError>;

    /// Restraint planes whose atoms deviate under `options`.
    ///
    /// # Errors
    ///
    /// As [`Self::nonplanar_aromatic_rings`].
    fn plane_restraint_outliers(
        &self,
        restraints: &[molframe_validate::PlaneRestraint],
        options: molframe_validate::PlanarityOptions,
    ) -> Result<molframe_validate::PlaneRestraintReport, molframe_validate::PlanarityError>;

    /// Ramachandran classification under `options`.
    ///
    /// # Errors
    ///
    /// Returns a Ramachandran error when the reference contours are unusable.
    fn ramachandran(
        &self,
        options: &molframe_validate::RamachandranOptions<'_>,
    ) -> Result<Vec<molframe_validate::RamachandranRecord>, molframe_validate::RamachandranError>;

    /// The subset of [`Self::ramachandran`] outside the allowed regions.
    ///
    /// # Errors
    ///
    /// As [`Self::ramachandran`].
    fn ramachandran_outliers(
        &self,
        options: &molframe_validate::RamachandranOptions<'_>,
    ) -> Result<Vec<molframe_validate::RamachandranRecord>, molframe_validate::RamachandranError>;

    /// Side-chain rotamers outside the `profile` under `options`.
    ///
    /// # Errors
    ///
    /// Returns a rotamer error for a residue the profile does not cover.
    fn rotamer_outliers(
        &self,
        provider: &dyn ComponentProvider,
        policy: &AnalysisPolicy,
        references: &molframe_validate::ReferenceLibrary,
        profile: &molframe_validate::RotamerProfile,
        options: molframe_validate::RotamerOptions,
    ) -> Result<molframe_validate::RotamerReport, molframe_validate::RotamerError>;

    /// Stereocentres whose configuration is ambiguous under `options`.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic when the chemistry cannot resolve a centre.
    fn chirality_outliers(
        &self,
        provider: &dyn ComponentProvider,
        policy: &AnalysisPolicy,
        options: molframe_validate::ChiralityOptions,
    ) -> Result<molframe_validate::ChiralityReport, molframe_core::diagnostic::Diagnostic>;

    /// Component atoms the dictionary defines but the structure omits.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic when the provider has no entry for a component.
    fn ccd_missing_atoms(
        &self,
        provider: &dyn ComponentProvider,
        policy: &AnalysisPolicy,
    ) -> Result<molframe_validate::CcdCompletenessReport, molframe_core::diagnostic::Diagnostic>;

    /// Nucleic-acid geometry under `policy`, with `roles` naming the atoms.
    ///
    /// # Errors
    ///
    /// Returns a geometry error for a nucleotide whose trace is incomplete.
    fn nucleic_acid_geometry(
        &self,
        provider: &dyn ComponentProvider,
        roles: &PolymerRoleProfile,
        policy: molframe_validate::NucleicGeometryPolicy,
    ) -> Result<
        Vec<molframe_validate::NucleicGeometryRecord>,
        molframe_validate::NucleicGeometryError,
    >;

    /// Per-chain completeness in `namespace`.
    ///
    /// # Errors
    ///
    /// Returns a completeness error when the structure has no chains.
    fn completeness(
        &self,
        namespace: Namespace,
    ) -> Result<Vec<molframe_validate::ChainCompleteness>, molframe_validate::CompletenessError>;

    /// Alternate-location occupancy sums in `namespace`.
    ///
    /// # Errors
    ///
    /// Returns an occupancy error for an inconsistent altloc partition.
    fn altloc_occupancy_sums(
        &self,
        namespace: Namespace,
        options: molframe_validate::AltlocOccupancyOptions,
    ) -> Result<molframe_validate::AltlocOccupancyReport, molframe_validate::AltlocOccupancyError>;

    /// Bond and angle deviations against the dictionary in `namespace`.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic when the provider has no reference geometry.
    fn reference_geometry(
        &self,
        provider: &dyn ComponentProvider,
        namespace: Namespace,
        options: molframe_validate::ReferenceGeometryOptions,
    ) -> Result<molframe_validate::ReferenceGeometryReport, molframe_core::diagnostic::Diagnostic>;

    /// The B-factor distribution over `selection`.
    ///
    /// # Errors
    ///
    /// Returns a B-factor error for an invalid deviation count or an empty
    /// selection.
    fn b_factor_distribution(
        &self,
        selection: &Selection,
        outlier_standard_deviations: f64,
    ) -> Result<molframe_validate::BFactorDistribution, molframe_validate::BFactorError>;

    /// TLS-model consistency across `groups`.
    ///
    /// # Errors
    ///
    /// Returns a B-factor error for a group that does not match the structure.
    fn tls_b_factor_consistency(
        &self,
        groups: &[molframe_validate::TlsGroup],
        maximum_absolute_deviation: f64,
        symmetry_tolerance: f64,
    ) -> Result<molframe_validate::TlsBFactorReport, molframe_validate::BFactorError>;
}

impl ValidationExt for Structure {
    fn clashes(
        &self,
        tolerance: f32,
        radius_set: RadiusSet,
    ) -> Result<molframe_validate::ClashTable, molframe_spatial::SpatialError> {
        self.clashes_with(
            tolerance,
            radius_set,
            SpatialBackend::Auto,
            &ExecutionContext::default(),
        )
    }

    fn clashes_with(
        &self,
        tolerance: f32,
        radius_set: RadiusSet,
        backend: SpatialBackend,
        context: &ExecutionContext,
    ) -> Result<molframe_validate::ClashTable, molframe_spatial::SpatialError> {
        molframe_validate::clashes(self.engine(), tolerance, radius_set, backend, context)
    }

    fn bond_length_deviations(&self, tolerance: f32) -> Vec<molframe_validate::BondDeviation> {
        molframe_validate::bond_length_deviations(self.engine(), tolerance)
    }

    fn ligand_geometry(&self, tolerance: f32) -> molframe_validate::LigandGeometryReport {
        molframe_validate::ligand_geometry(self.engine(), tolerance)
    }

    fn ligand_geometry_outliers(&self, tolerance: f32) -> Vec<molframe_validate::BondDeviation> {
        molframe_validate::ligand_geometry_outliers(self.engine(), tolerance)
    }

    fn quality_flags(&self) -> Vec<molframe_validate::QualityFlag> {
        molframe_validate::quality_flags(self.engine())
    }

    fn overvalent_atoms(&self) -> Vec<molframe_validate::ValenceError> {
        molframe_validate::overvalent_atoms(self.engine())
    }

    fn cis_peptides(
        &self,
        threshold_degrees: f64,
    ) -> Result<Vec<molframe_validate::CisPeptide>, molframe_core::diagnostic::Diagnostic> {
        molframe_validate::cis_peptides(self.engine(), threshold_degrees)
    }

    fn nonplanar_aromatic_rings(
        &self,
        options: molframe_validate::PlanarityOptions,
    ) -> Result<Vec<molframe_validate::PlanarityFlag>, molframe_validate::PlanarityError> {
        molframe_validate::nonplanar_aromatic_rings(self.engine(), options)
    }

    fn plane_restraint_outliers(
        &self,
        restraints: &[molframe_validate::PlaneRestraint],
        options: molframe_validate::PlanarityOptions,
    ) -> Result<molframe_validate::PlaneRestraintReport, molframe_validate::PlanarityError> {
        molframe_validate::plane_restraint_outliers(self.engine(), restraints, options)
    }

    fn ramachandran(
        &self,
        options: &molframe_validate::RamachandranOptions<'_>,
    ) -> Result<Vec<molframe_validate::RamachandranRecord>, molframe_validate::RamachandranError>
    {
        molframe_validate::ramachandran(self.engine(), options)
    }

    fn ramachandran_outliers(
        &self,
        options: &molframe_validate::RamachandranOptions<'_>,
    ) -> Result<Vec<molframe_validate::RamachandranRecord>, molframe_validate::RamachandranError>
    {
        molframe_validate::ramachandran_outliers(self.engine(), options)
    }

    fn rotamer_outliers(
        &self,
        provider: &dyn ComponentProvider,
        policy: &AnalysisPolicy,
        references: &molframe_validate::ReferenceLibrary,
        profile: &molframe_validate::RotamerProfile,
        options: molframe_validate::RotamerOptions,
    ) -> Result<molframe_validate::RotamerReport, molframe_validate::RotamerError> {
        molframe_validate::rotamer_outliers(
            self.engine(),
            provider,
            policy,
            references,
            profile,
            options,
        )
    }

    fn chirality_outliers(
        &self,
        provider: &dyn ComponentProvider,
        policy: &AnalysisPolicy,
        options: molframe_validate::ChiralityOptions,
    ) -> Result<molframe_validate::ChiralityReport, molframe_core::diagnostic::Diagnostic> {
        molframe_validate::chirality_outliers(self.engine(), provider, policy, options)
    }

    fn ccd_missing_atoms(
        &self,
        provider: &dyn ComponentProvider,
        policy: &AnalysisPolicy,
    ) -> Result<molframe_validate::CcdCompletenessReport, molframe_core::diagnostic::Diagnostic>
    {
        molframe_validate::ccd_missing_atoms(self.engine(), provider, policy)
    }

    fn nucleic_acid_geometry(
        &self,
        provider: &dyn ComponentProvider,
        roles: &PolymerRoleProfile,
        policy: molframe_validate::NucleicGeometryPolicy,
    ) -> Result<
        Vec<molframe_validate::NucleicGeometryRecord>,
        molframe_validate::NucleicGeometryError,
    > {
        molframe_validate::nucleic_acid_geometry(self.engine(), provider, roles, policy)
    }

    fn completeness(
        &self,
        namespace: Namespace,
    ) -> Result<Vec<molframe_validate::ChainCompleteness>, molframe_validate::CompletenessError>
    {
        molframe_validate::completeness(self.engine(), namespace)
    }

    fn altloc_occupancy_sums(
        &self,
        namespace: Namespace,
        options: molframe_validate::AltlocOccupancyOptions,
    ) -> Result<molframe_validate::AltlocOccupancyReport, molframe_validate::AltlocOccupancyError>
    {
        molframe_validate::altloc_occupancy_sums(self.engine(), namespace, options)
    }

    fn reference_geometry(
        &self,
        provider: &dyn ComponentProvider,
        namespace: Namespace,
        options: molframe_validate::ReferenceGeometryOptions,
    ) -> Result<molframe_validate::ReferenceGeometryReport, molframe_core::diagnostic::Diagnostic>
    {
        molframe_validate::reference_geometry(self.engine(), provider, namespace, options)
    }

    fn b_factor_distribution(
        &self,
        selection: &Selection,
        outlier_standard_deviations: f64,
    ) -> Result<molframe_validate::BFactorDistribution, molframe_validate::BFactorError> {
        molframe_validate::b_factor_distribution(
            self.engine(),
            selection.atom_selection(),
            outlier_standard_deviations,
        )
    }

    fn tls_b_factor_consistency(
        &self,
        groups: &[molframe_validate::TlsGroup],
        maximum_absolute_deviation: f64,
        symmetry_tolerance: f64,
    ) -> Result<molframe_validate::TlsBFactorReport, molframe_validate::BFactorError> {
        molframe_validate::tls_b_factor_consistency(
            self.engine(),
            groups,
            maximum_absolute_deviation,
            symmetry_tolerance,
        )
    }
}
