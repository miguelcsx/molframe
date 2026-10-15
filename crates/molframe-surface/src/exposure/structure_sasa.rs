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
    AppliedRadiusPolicy, HeteroPolicy, MetalFallback, MetalTables, PolicyWordError, RadiusPolicy,
    RadiusSource, TableError, TableRequest, TableStore, VdwRadii, policy_radii,
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
/// There is no `Default`: the radii are a table, and a table is chosen, not
/// assumed. [`SasaOptions::resolve`] takes the choice and fills the rest with
/// Shrake-Rupley at 960 points, a 1.4 Å probe, the polymer alone (no waters,
/// ligands or ions), the vdW metal fallback and the first model.
/// [`SasaOptions::applied_policy`] reports which table, version, source and
/// licence a calculation used.
#[derive(Clone, Debug, PartialEq)]
pub struct SasaOptions {
    /// The sampling method.
    pub method: SasaMethod,
    /// Probe radius in ångström.
    pub probe: f32,
    /// The radius table in force.
    pub radii: RadiusSource,
    /// The tables behind the metal fallback; absent when `metals` is `None`.
    pub metal_tables: Option<MetalTables>,
    /// Which heteroatom classes are measured and occlude.
    pub hetero: HeteroPolicy,
    /// Where a metal the radius set has no value for takes its radius.
    pub metals: MetalFallback,
    /// Which model to measure.
    pub model: ModelIndex,
}

impl SasaOptions {
    /// The default method, probe, heteroatom and metal choices over the
    /// radius table `radii` selects.
    ///
    /// # Errors
    ///
    /// [`TableError`] naming how to supply whichever table is missing.
    pub fn resolve(radii: &TableRequest<VdwRadii>, store: &TableStore) -> Result<Self, TableError> {
        Self::assemble(
            SasaMethod::default(),
            DEFAULT_PROBE,
            HeteroPolicy::POLYMER,
            MetalFallback::Vdw,
            radii,
            store,
        )
    }

    fn assemble(
        method: SasaMethod,
        probe: f32,
        hetero: HeteroPolicy,
        metals: MetalFallback,
        radii: &TableRequest<VdwRadii>,
        store: &TableStore,
    ) -> Result<Self, TableError> {
        let policy = RadiusPolicy::resolve(radii, metals, hetero, store)?;
        Ok(Self {
            method,
            probe,
            radii: policy.source,
            metal_tables: policy.metal_tables,
            hetero,
            metals,
            model: ModelIndex::new(0),
        })
    }

    /// The radius assignment these options describe.
    #[must_use]
    pub fn radius_policy(&self) -> RadiusPolicy {
        RadiusPolicy {
            source: self.radii.clone(),
            metals: self.metals,
            metal_tables: self.metal_tables.clone(),
            hetero: self.hetero,
        }
    }

    /// The radius policy with every table it resolves to named, as a result
    /// should record it.
    #[must_use]
    pub fn applied_policy(&self) -> AppliedRadiusPolicy {
        self.radius_policy().applied()
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
    /// A radius table could not be supplied.
    #[error(transparent)]
    Table(#[from] TableError),
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
    /// and [`MetalFallback`] parsing; `radii` selects the radius table.
    ///
    /// # Errors
    ///
    /// [`SasaWordError`] naming the word that was not recognised, or the table
    /// that could not be supplied.
    pub fn from_words(
        method: Option<&str>,
        samples: Option<u16>,
        probe: Option<f32>,
        hetero: Option<&str>,
        metals: Option<&str>,
        radii: &TableRequest<VdwRadii>,
        store: &TableStore,
    ) -> Result<Self, SasaWordError> {
        let method = SasaMethod::named(
            match method {
                Some(name) => name,
                None => "shrake-rupley",
            },
            samples,
        )?;
        let hetero = match hetero {
            Some(words) => words.parse()?,
            None => HeteroPolicy::POLYMER,
        };
        let metals = match metals {
            Some(words) => words.parse()?,
            None => MetalFallback::Vdw,
        };
        let probe = match probe {
            Some(value) => value,
            None => DEFAULT_PROBE,
        };
        Ok(Self::assemble(method, probe, hetero, metals, radii, store)?)
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
