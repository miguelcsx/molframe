//! Structure-first solvent-accessible surface area, as methods.
//!
//! Unlike the other extension traits, this one is not a pure forwarder: the
//! surface kernels take positions and radii, so the method decides what to
//! gather and what to do about an atom whose element has no radius in the
//! chosen set. The decision is stated here once. Such an atom gets a `NaN`
//! area and, because it has no radius, does not occlude its neighbours.

use crate::ExecutionContext;
use crate::structure::Structure;
use molframe_chem::{RadiusSet, atom_radii};
use molframe_core::diagnostic::Diagnostic;
use molframe_core::index::ModelIndex;
use molframe_surface::SasaError;

/// The sampling method behind a SASA calculation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SasaMethod {
    /// Shrake-Rupley with this many test points per atom.
    ShrakeRupley {
        /// Test points per atom sphere.
        points: u16,
    },
    /// Lee-Richards with this many slices per atom sphere.
    LeeRichards {
        /// Slices per atom sphere.
        slices: u16,
    },
}

/// Options for [`SurfaceExt::sasa_with`].
///
/// The default is Shrake-Rupley with 960 points, a 1.4 Å probe, Bondi radii
/// and the first model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SasaOptions {
    /// The sampling method.
    pub method: SasaMethod,
    /// Probe radius in ångström.
    pub probe: f32,
    /// Which radius table to use.
    pub radii: RadiusSet,
    /// Which model to measure.
    pub model: ModelIndex,
}

impl Default for SasaOptions {
    fn default() -> Self {
        Self {
            method: SasaMethod::ShrakeRupley { points: 960 },
            probe: 1.4,
            radii: RadiusSet::Bondi,
            model: ModelIndex::new(0),
        }
    }
}

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
    /// the chosen set get `NaN`.
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
        let radii = atom_radii(engine, options.radii);
        // An atom without a radius cannot occlude and has no area: give it a
        // zero-radius sphere to measure, then report NaN for it.
        let measured: Vec<f32> = radii
            .iter()
            .map(|radius| if radius.is_nan() { 0.0 } else { *radius })
            .collect();
        let mut areas = match options.method {
            SasaMethod::ShrakeRupley { points } => molframe_surface::shrake_rupley(
                positions,
                &measured,
                options.probe,
                points,
                context,
            )?,
            SasaMethod::LeeRichards { slices } => molframe_surface::lee_richards(
                positions,
                &measured,
                options.probe,
                slices,
                context,
            )?,
        };
        for (area, radius) in areas.iter_mut().zip(&radii) {
            if radius.is_nan() {
                *area = f64::NAN;
            }
        }
        Ok(areas)
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
