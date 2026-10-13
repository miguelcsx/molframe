//! One declarative option set for solvent-accessible area of a structure.
//!
//! The kernels in this crate take positions and radii. A structure-level
//! caller has to decide which atoms get a radius, and every front end (the
//! facade, the analysis crate, the command line and Python) must decide it the
//! same way, so the decision is stated once, here.
//!
//! The default measures the polymer alone, as NACCESS and `FreeSASA` do: waters,
//! ligands and ions get `NaN` and do not occlude unless [`SasaOptions::hetero`]
//! brings them in. [`SasaOptions::metals`] says where a metal the radius set
//! lacks takes its radius from. An atom without a positive finite radius is not
//! an obstacle: even a zero-radius sphere grows by the probe and would shield
//! its neighbours, so such atoms are left out of the kernel entirely and their
//! area is `NaN`. Both sampling methods read the one radius vector, so
//! Shrake-Rupley and Lee-Richards always measure the same set of atoms.

use crate::{SasaError, lee_richards, shrake_rupley};
use molframe_chem::{
    AppliedRadiusPolicy, CustomRadii, HeteroPolicy, MetalFallback, PolicyWordError, RadiusPolicy,
    RadiusSet, RadiusSource, RadiusTableError, policy_radii,
};
use molframe_core::ExecutionContext;
use molframe_core::index::ModelIndex;
use molframe_core::structure::Structure;

/// Test points per atom of the default Shrake-Rupley sampling.
pub const DEFAULT_POINTS: u16 = 960;
/// Slices per atom of the default Lee-Richards sampling.
pub const DEFAULT_SLICES: u16 = 20;
/// The default probe radius in ångström.
pub const DEFAULT_PROBE: f32 = 1.4;

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

impl Default for SasaMethod {
    fn default() -> Self {
        Self::ShrakeRupley {
            points: DEFAULT_POINTS,
        }
    }
}

/// Options for a structure-level SASA calculation.
///
/// The default is Shrake-Rupley with 960 points, a 1.4 Å probe, Bondi radii,
/// the polymer alone (no waters, ligands or ions), the vdW metal fallback and
/// the first model.
///
/// The radius source is an argument, not a constant: [`SasaOptions::radii`]
/// holds either one of the built-in sets or a table the caller supplied (see
/// [`SasaOptions::with_radii_table`]), and [`SasaOptions::applied_policy`]
/// reports which table, version, source and licence a calculation used.
#[derive(Clone, Debug, PartialEq)]
pub struct SasaOptions {
    /// The sampling method.
    pub method: SasaMethod,
    /// Probe radius in ångström.
    pub probe: f32,
    /// Which radius table to use: a built-in set or the caller's own.
    pub radii: RadiusSource,
    /// Which heteroatom classes are measured and occlude.
    pub hetero: HeteroPolicy,
    /// Where a metal the radius set has no value for takes its radius.
    pub metals: MetalFallback,
    /// Which model to measure.
    pub model: ModelIndex,
}

impl Default for SasaOptions {
    fn default() -> Self {
        Self {
            method: SasaMethod::default(),
            probe: DEFAULT_PROBE,
            radii: RadiusSource::Set(RadiusSet::Bondi),
            hetero: HeteroPolicy::POLYMER,
            metals: MetalFallback::Vdw,
            model: ModelIndex::new(0),
        }
    }
}

impl SasaOptions {
    /// The radius assignment these options describe.
    #[must_use]
    pub fn radius_policy(&self) -> RadiusPolicy {
        RadiusPolicy {
            source: self.radii.clone(),
            metals: self.metals,
            hetero: self.hetero,
        }
    }

    /// The radius policy with every table it resolves to named, as a result
    /// should record it.
    #[must_use]
    pub fn applied_policy(&self) -> AppliedRadiusPolicy {
        self.radius_policy().applied()
    }

    /// These options with the radii taken from a caller's table, given as
    /// text in the format of [`CustomRadii::parse`].
    ///
    /// # Errors
    ///
    /// [`RadiusTableError`] for a table that does not validate.
    pub fn with_radii_table(self, text: &str) -> Result<Self, RadiusTableError> {
        Ok(Self {
            radii: RadiusSource::from(CustomRadii::parse(text)?),
            ..self
        })
    }

    /// The radius of every atom of `structure`, indexed by atom. `NaN` marks
    /// an atom that is left out: excluded by `hetero`, or without a radius.
    #[must_use]
    pub fn atom_radii(&self, structure: &Structure) -> Vec<f32> {
        policy_radii(structure, &self.radius_policy())
    }
}

/// A word naming a method, heteroatom class or metal fallback that does not exist.
#[derive(Debug, thiserror::Error)]
pub enum SasaWordError {
    /// The method name is neither `shrake-rupley` nor `lee-richards`.
    #[error("unknown SASA method `{0}`; expected shrake-rupley or lee-richards")]
    Method(String),
    /// A heteroatom or metal word was not recognised.
    #[error(transparent)]
    Policy(#[from] PolicyWordError),
    /// The sample count does not fit the kernel's range.
    #[error("the sample count must be between 1 and 65535")]
    Samples,
}

impl SasaMethod {
    /// The method called `name` (`shrake-rupley` or `lee-richards`), with
    /// `samples` points or slices, or the method's own default when `None`.
    ///
    /// # Errors
    ///
    /// [`SasaWordError::Method`] for another name, [`SasaWordError::Samples`]
    /// for a count of zero.
    pub fn named(name: &str, samples: Option<u16>) -> Result<Self, SasaWordError> {
        if samples == Some(0) {
            return Err(SasaWordError::Samples);
        }
        match name {
            "shrake-rupley" | "sr" => Ok(Self::ShrakeRupley {
                points: match samples {
                    Some(count) => count,
                    None => DEFAULT_POINTS,
                },
            }),
            "lee-richards" | "lr" => Ok(Self::LeeRichards {
                slices: match samples {
                    Some(count) => count,
                    None => DEFAULT_SLICES,
                },
            }),
            other => Err(SasaWordError::Method(other.to_owned())),
        }
    }
}

impl SasaOptions {
    /// Options built from the words a front end receives; every `None` keeps
    /// the default. `hetero` and `metals` use the spellings of [`HeteroPolicy`]
    /// and [`MetalFallback`] parsing.
    ///
    /// # Errors
    ///
    /// [`SasaWordError`] naming the word that was not recognised.
    pub fn from_words(
        method: Option<&str>,
        samples: Option<u16>,
        probe: Option<f32>,
        hetero: Option<&str>,
        metals: Option<&str>,
    ) -> Result<Self, SasaWordError> {
        let defaults = Self::default();
        Ok(Self {
            method: match method {
                Some(name) => SasaMethod::named(name, samples)?,
                None => SasaMethod::named("shrake-rupley", samples)?,
            },
            probe: match probe {
                Some(value) => value,
                None => defaults.probe,
            },
            hetero: match hetero {
                Some(words) => words.parse()?,
                None => defaults.hetero,
            },
            metals: match metals {
                Some(words) => words.parse()?,
                None => defaults.metals,
            },
            ..defaults
        })
    }
}

/// Per-atom area for `positions` with one radius per atom.
///
/// An entry of `radii` that is not finite and positive removes that atom from
/// the calculation and yields `NaN` for it.
///
/// # Errors
///
/// [`SasaError`] for an invalid probe or sampling density, or a refused budget.
pub fn measure_atoms(
    positions: &[[f32; 3]],
    radii: &[f32],
    options: &SasaOptions,
    context: &ExecutionContext,
) -> Result<Vec<f64>, SasaError> {
    let kept: Vec<usize> = (0..radii.len().min(positions.len()))
        .filter(|&atom| radii[atom].is_finite() && radii[atom] > 0.0)
        .collect();
    let points: Vec<[f32; 3]> = kept.iter().map(|&atom| positions[atom]).collect();
    let spheres: Vec<f32> = kept.iter().map(|&atom| radii[atom]).collect();
    let measured = match options.method {
        SasaMethod::ShrakeRupley { points: samples } => {
            shrake_rupley(&points, &spheres, options.probe, samples, context)?
        }
        SasaMethod::LeeRichards { slices } => {
            lee_richards(&points, &spheres, options.probe, slices, context)?
        }
    };
    let mut areas = vec![f64::NAN; radii.len()];
    for (&atom, area) in kept.iter().zip(measured) {
        areas[atom] = area;
    }
    Ok(areas)
}

/// Per-atom areas together with the radius policy that produced them.
#[derive(Clone, Debug, PartialEq)]
pub struct SasaReport {
    /// Area of every atom, indexed by atom; `NaN` where the atom is left out.
    pub areas: Vec<f64>,
    /// The tables and choices the radii came from.
    pub policy: AppliedRadiusPolicy,
}

/// Like [`structure_atom_areas`], and states which radius policy applied.
///
/// # Errors
///
/// See [`measure_atoms`].
pub fn structure_sasa_report(
    structure: &Structure,
    positions: &[[f32; 3]],
    options: &SasaOptions,
    context: &ExecutionContext,
) -> Result<SasaReport, SasaError> {
    Ok(SasaReport {
        areas: structure_atom_areas(structure, positions, options, context)?,
        policy: options.applied_policy(),
    })
}

/// Per-atom area of `structure` at `positions`, indexed by atom.
///
/// # Errors
///
/// See [`measure_atoms`].
pub fn structure_atom_areas(
    structure: &Structure,
    positions: &[[f32; 3]],
    options: &SasaOptions,
    context: &ExecutionContext,
) -> Result<Vec<f64>, SasaError> {
    measure_atoms(positions, &options.atom_radii(structure), options, context)
}

#[cfg(test)]
#[path = "structure_sasa_tests.rs"]
mod tests;
