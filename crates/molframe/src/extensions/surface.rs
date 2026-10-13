//! Structure-first solvent-accessible surface area, as methods.
//!
//! Unlike the other extension traits, this one is not a pure forwarder: the
//! surface kernels take positions and radii, so the method gathers the model's
//! coordinates and hands them to [`molframe_surface::structure_atom_areas`].
//! That function owns the decisions every front end shares: the polymer alone
//! is measured by default, waters, ligands and ions come in through
//! [`SasaOptions::hetero`], and an atom with no positive radius gets a `NaN`
//! area and does not occlude its neighbours.

use crate::ExecutionContext;
use crate::structure::Structure;
use molframe_core::diagnostic::Diagnostic;
use molframe_core::index::ModelIndex;
use molframe_surface::SasaError;
pub use molframe_surface::{SasaMethod, SasaOptions, SasaReport};

/// Why a structure-level SASA call failed.
#[derive(Debug, thiserror::Error)]
pub enum SurfaceExtError {
    /// The structure is a ragged ensemble or has no such model. The message
    /// names `Structure::ragged_models`, which is the way forward.
    #[error("{0}")]
    Structure(Box<Diagnostic>),
    /// The surface kernel refused the input.
    #[error(transparent)]
    Sasa(#[from] SasaError),
}

/// Surface-area operations over a structure.
///
/// # Examples
///
/// ```
/// # #[cfg(all(feature = "analysis", feature = "chemistry"))]
/// # {
/// use molframe::prelude::*;
///
/// # const PDB: &str = "\
/// # ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N
/// # ATOM      2  CA  ALA A   1      12.560   6.195  -6.504  1.00  0.00           C
/// # END
/// # ";
/// let (structure, _) = read_bytes(PDB.into(), Some("1abc.pdb"), &ReadOptions::new())?;
/// let areas = structure.sasa()?;
/// assert_eq!(areas.len(), structure.atom_count() as usize);
/// let per_residue = structure.residue_sasa()?;
/// assert_eq!(per_residue.len(), 1);
///
/// // The README quick start: a selection, an explicit memory ceiling, contacts.
/// let alpha = structure.select("name CA", &AnalysisPolicy::default())?;
/// assert_eq!(alpha.to_coordinates(ModelIndex::new(0)).len(), 1);
/// let context = ExecutionContext::builder()
///     .memory_budget(MemoryBudget::new(512 << 20)?)
///     .build()?;
/// let contacts = structure.atom_contacts_with(4.5, SpatialBackend::Auto, &context)?;
/// assert!(!contacts.is_empty());
/// # }
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub trait SurfaceExt {
    /// Per-atom SASA in square ångström with the default options.
    ///
    /// # Errors
    ///
    /// See [`Self::sasa_with`].
    fn sasa(&self) -> Result<Vec<f64>, SurfaceExtError>;

    /// Per-atom SASA, indexed by atom. Atoms whose element has no radius in
    /// the chosen set, and heteroatoms the options leave out, get `NaN`.
    ///
    /// # Errors
    ///
    /// [`SurfaceExtError::Structure`] for a ragged ensemble or a missing model;
    /// [`SurfaceExtError::Sasa`] for an invalid probe or sampling density.
    fn sasa_with(
        &self,
        options: &SasaOptions,
        context: &ExecutionContext,
    ) -> Result<Vec<f64>, SurfaceExtError>;

    /// Per-atom SASA together with the radius policy that applied: the radius
    /// set or caller-supplied table, the metal tables, and the hetero choice,
    /// each by name, version, source and licence.
    ///
    /// # Errors
    ///
    /// See [`Self::sasa_with`].
    fn sasa_report_with(
        &self,
        options: &SasaOptions,
        context: &ExecutionContext,
    ) -> Result<SasaReport, SurfaceExtError>;

    /// Per-residue SASA with the default options, indexed by residue.
    ///
    /// # Errors
    ///
    /// See [`Self::sasa_with`].
    fn residue_sasa(&self) -> Result<Vec<f64>, SurfaceExtError>;

    /// Per-residue SASA: the sum of the finite atom areas of each residue.
    ///
    /// # Errors
    ///
    /// See [`Self::sasa_with`].
    fn residue_sasa_with(
        &self,
        options: &SasaOptions,
        context: &ExecutionContext,
    ) -> Result<Vec<f64>, SurfaceExtError>;
}

impl SurfaceExt for Structure {
    fn sasa(&self) -> Result<Vec<f64>, SurfaceExtError> {
        self.sasa_with(&SasaOptions::default(), &ExecutionContext::default())
    }

    fn sasa_with(
        &self,
        options: &SasaOptions,
        context: &ExecutionContext,
    ) -> Result<Vec<f64>, SurfaceExtError> {
        let engine = self.engine();
        engine
            .require_resolved_model()
            .map_err(|finding| SurfaceExtError::Structure(Box::new(finding)))?;
        let positions = self
            .model_coordinates(options.model)
            .ok_or_else(|| missing_model(self, options.model))?;
        Ok(molframe_surface::structure_atom_areas(
            engine, positions, options, context,
        )?)
    }

    fn sasa_report_with(
        &self,
        options: &SasaOptions,
        context: &ExecutionContext,
    ) -> Result<SasaReport, SurfaceExtError> {
        Ok(SasaReport {
            areas: self.sasa_with(options, context)?,
            policy: options.applied_policy(),
        })
    }

    fn residue_sasa(&self) -> Result<Vec<f64>, SurfaceExtError> {
        self.residue_sasa_with(&SasaOptions::default(), &ExecutionContext::default())
    }

    fn residue_sasa_with(
        &self,
        options: &SasaOptions,
        context: &ExecutionContext,
    ) -> Result<Vec<f64>, SurfaceExtError> {
        let areas = self.sasa_with(options, context)?;
        let mut totals = vec![0.0_f64; self.residues().len()];
        for atom in self.atoms() {
            let area = match areas.get(atom.index().as_usize()).copied() {
                Some(area) => area,
                None => f64::NAN,
            };
            if !area.is_finite() {
                continue;
            }
            if let Some(slot) = atom
                .residue()
                .and_then(|residue| totals.get_mut(residue.index().as_usize()))
            {
                *slot += area;
            }
        }
        Ok(totals)
    }
}

fn missing_model(structure: &Structure, model: ModelIndex) -> SurfaceExtError {
    // `ensemble_positions` is the one place that words a missing model, so the
    // façade reuses its diagnostic rather than inventing a second message.
    match structure.engine().ensemble_positions(Some(&[model])) {
        Err(finding) => SurfaceExtError::Structure(Box::new(finding)),
        Ok(_) => SurfaceExtError::Sasa(SasaError::NoPoints),
    }
}
