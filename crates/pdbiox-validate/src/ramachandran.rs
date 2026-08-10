//! Backbone conformation classified by versioned empirical φ/ψ grids.
//!
//! Coordinates are traversed once, while each measured torsion pair is assessed
//! against the explicitly selected basin grids. For `r` residues and `b`
//! configured basins this costs `O(r * b * log n)`, where `n` is the larger grid
//! axis. No built-in residue names, rectangular basins, thresholds or reference
//! releases participate in classification.

use crate::distributions::{ReferenceAssessment, ReferenceError, ReferenceLibrary};
use pdbiox_chem::PolymerAtomRole;
use pdbiox_core::Diagnostic;
use pdbiox_core::index::ResidueIndex;
use pdbiox_core::structure::{ResidueRef, Structure};
use std::collections::BTreeSet;

/// The named conformational basin containing a φ/ψ observation.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum RamachandranRegion {
    /// The right-handed α-helix basin.
    AlphaHelixRight,
    /// The extended β-strand basin.
    BetaSheet,
    /// The left-handed α-helix basin.
    AlphaHelixLeft,
    /// Empirical support does not reach the configured minimum.
    Outlier,
}

/// Associates a conformational basin with one grid in a reference library.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RamachandranBasin {
    region: RamachandranRegion,
    distribution: Box<str>,
}

impl RamachandranBasin {
    /// Selects a named grid for a conformational basin.
    ///
    /// # Errors
    ///
    /// Refuses an empty distribution name or the outlier pseudo-region.
    pub fn new(
        region: RamachandranRegion,
        distribution: impl Into<Box<str>>,
    ) -> Result<Self, RamachandranError> {
        if region == RamachandranRegion::Outlier {
            return Err(RamachandranError::OutlierBasin);
        }
        let distribution = distribution.into();
        if distribution.trim().is_empty() {
            return Err(RamachandranError::EmptyDistribution);
        }
        Ok(Self {
            region,
            distribution,
        })
    }

    /// Basin assigned to observations supported by this grid.
    #[must_use]
    pub const fn region(&self) -> RamachandranRegion {
        self.region
    }

    /// Grid name inside the selected reference library.
    #[must_use]
    pub fn distribution(&self) -> &str {
        &self.distribution
    }
}

/// Explicit data and decision policy for Ramachandran classification.
#[derive(Clone, Debug, PartialEq)]
pub struct RamachandranOptions<'a> {
    references: &'a ReferenceLibrary,
    basins: Vec<RamachandranBasin>,
    minimum_probability: f64,
}

impl<'a> RamachandranOptions<'a> {
    /// Creates a reproducible classification policy.
    ///
    /// The highest-probability configured basin wins. Input order resolves an
    /// exact tie, so callers can make tie-breaking stable and visible.
    ///
    /// # Errors
    ///
    /// Refuses no basins, duplicate regions, or a non-finite probability
    /// outside the closed unit interval.
    pub fn new(
        references: &'a ReferenceLibrary,
        basins: impl IntoIterator<Item = RamachandranBasin>,
        minimum_probability: f64,
    ) -> Result<Self, RamachandranError> {
        if !minimum_probability.is_finite() || !(0.0..=1.0).contains(&minimum_probability) {
            return Err(RamachandranError::Probability);
        }
        let basins: Vec<_> = basins.into_iter().collect();
        if basins.is_empty() {
            return Err(RamachandranError::NoBasins);
        }
        let unique: BTreeSet<_> = basins.iter().map(RamachandranBasin::region).collect();
        if unique.len() != basins.len() {
            return Err(RamachandranError::DuplicateBasin);
        }
        let distributions: BTreeSet<_> =
            basins.iter().map(RamachandranBasin::distribution).collect();
        if distributions.len() != basins.len() {
            return Err(RamachandranError::DuplicateDistribution);
        }
        for basin in &basins {
            references.validate_grid(basin.distribution())?;
        }
        Ok(Self {
            references,
            basins,
            minimum_probability,
        })
    }

    /// Exact reference library used by the classifier.
    #[must_use]
    pub const fn references(&self) -> &ReferenceLibrary {
        self.references
    }

    /// Configured basin grids in deterministic tie-breaking order.
    #[must_use]
    pub fn basins(&self) -> &[RamachandranBasin] {
        &self.basins
    }

    /// Minimum normalized grid-bin mass accepted as supported.
    #[must_use]
    pub const fn minimum_probability(&self) -> f64 {
        self.minimum_probability
    }
}

/// A residue placed on a versioned Ramachandran reference grid.
#[derive(Clone, Debug, PartialEq)]
pub struct RamachandranRecord {
    /// The residue described.
    pub residue: ResidueIndex,
    /// φ in degrees.
    pub phi: f64,
    /// ψ in degrees.
    pub psi: f64,
    /// Highest-support basin, or outlier below the configured threshold.
    pub region: RamachandranRegion,
    /// Reference identity, release, grid name and normalized empirical support.
    pub reference: ReferenceAssessment,
}

/// Invalid policy, reference grid or observation.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum RamachandranError {
    /// Assessable residue counts exceed public coverage counters.
    #[error("Ramachandran residue count exceeds u32 coverage limits")]
    CoverageOverflow,
    /// No basin grid was configured.
    #[error("at least one Ramachandran basin grid is required")]
    NoBasins,
    /// A basin occurs more than once in the policy.
    #[error("Ramachandran basin regions must be unique")]
    DuplicateBasin,
    /// One grid was assigned to more than one basin.
    #[error("Ramachandran basin distributions must be unique")]
    DuplicateDistribution,
    /// The outlier pseudo-region cannot own a reference grid.
    #[error("the outlier region cannot be configured as a basin")]
    OutlierBasin,
    /// A selected distribution name is empty.
    #[error("Ramachandran distribution names must not be empty")]
    EmptyDistribution,
    /// Minimum probability is not finite or lies outside zero to one.
    #[error("minimum Ramachandran probability must be finite and between zero and one")]
    Probability,
    /// The selected reference data cannot assess the observation.
    #[error(transparent)]
    Reference(#[from] ReferenceError),
    /// Polymer atom roles are absent or ambiguous.
    #[error("polymer atom-role projection failed: {0}")]
    Roles(Diagnostic),
}

/// Classifies a finite φ/ψ pair against explicitly selected empirical grids.
///
/// # Errors
///
/// Returns an error when a selected distribution is absent, is not a grid, or
/// either observation is non-finite.
pub fn classify(
    phi: f64,
    psi: f64,
    options: &RamachandranOptions<'_>,
) -> Result<(RamachandranRegion, ReferenceAssessment), RamachandranError> {
    let mut best: Option<(RamachandranRegion, ReferenceAssessment)> = None;
    for basin in options.basins() {
        let assessment = options
            .references()
            .assess_pair(basin.distribution(), phi, psi)?;
        let replace = match &best {
            Some((_, current)) => assessment.probability > current.probability,
            None => true,
        };
        if replace {
            best = Some((basin.region(), assessment));
        }
    }
    let Some((region, assessment)) = best else {
        return Err(RamachandranError::NoBasins);
    };
    let region = if assessment.probability >= options.minimum_probability() {
        region
    } else {
        RamachandranRegion::Outlier
    };
    Ok((region, assessment))
}

/// Places every residue with defined φ and ψ on explicit reference grids.
///
/// Chain ends and residues across breaks are omitted. Records are ordered by
/// residue index.
///
/// # Errors
///
/// Returns the first reference-data error encountered.
pub fn ramachandran(
    structure: &Structure,
    options: &RamachandranOptions<'_>,
) -> Result<Vec<RamachandranRecord>, RamachandranError> {
    crate::backbone::require_polymer_roles(structure).map_err(RamachandranError::Roles)?;
    let mut records = Vec::new();
    for chain in structure.data().chains() {
        let residues: Vec<ResidueRef<'_>> = chain.residues().collect();
        for position in 0..residues.len() {
            let previous = position.checked_sub(1).and_then(|p| residues.get(p));
            let current = residues[position];
            let next = residues.get(position + 1);
            let (Some(previous), Some(next)) = (previous, next) else {
                continue;
            };
            let (Some(phi), Some(psi)) = (
                phi(structure, *previous, current)?,
                psi(structure, current, *next)?,
            ) else {
                continue;
            };
            let phi = pdbiox_geom::degrees(phi);
            let psi = pdbiox_geom::degrees(psi);
            let (region, reference) = classify(phi, psi, options)?;
            records.push(RamachandranRecord {
                residue: current.index(),
                phi,
                psi,
                region,
                reference,
            });
        }
    }
    records.sort_by_key(|record| record.residue.get());
    Ok(records)
}

/// Returns only observations below the configured empirical-support threshold.
///
/// # Errors
///
/// Returns the first reference-data error encountered.
pub fn ramachandran_outliers(
    structure: &Structure,
    options: &RamachandranOptions<'_>,
) -> Result<Vec<RamachandranRecord>, RamachandranError> {
    Ok(ramachandran(structure, options)?
        .into_iter()
        .filter(|record| record.region == RamachandranRegion::Outlier)
        .collect())
}

fn phi(
    structure: &Structure,
    previous: ResidueRef<'_>,
    current: ResidueRef<'_>,
) -> Result<Option<f64>, RamachandranError> {
    let previous_carbon = crate::backbone::role_atom(
        structure,
        previous,
        PolymerAtomRole::PROTEIN_CARBONYL_CARBON,
    )
    .map_err(RamachandranError::Roles)?;
    let nitrogen =
        crate::backbone::role_atom(structure, current, PolymerAtomRole::PROTEIN_NITROGEN)
            .map_err(RamachandranError::Roles)?;
    let alpha =
        crate::backbone::role_atom(structure, current, PolymerAtomRole::PROTEIN_ALPHA_CARBON)
            .map_err(RamachandranError::Roles)?;
    let carbon =
        crate::backbone::role_atom(structure, current, PolymerAtomRole::PROTEIN_CARBONYL_CARBON)
            .map_err(RamachandranError::Roles)?;
    let (Some(previous_carbon), Some(nitrogen), Some(alpha), Some(carbon)) =
        (previous_carbon, nitrogen, alpha, carbon)
    else {
        return Ok(None);
    };
    if !crate::backbone::peptide_bonded(structure, previous_carbon, nitrogen) {
        return Ok(None);
    }
    let (Some(previous_carbon), Some(nitrogen), Some(alpha), Some(carbon)) = (
        previous_carbon.position(),
        nitrogen.position(),
        alpha.position(),
        carbon.position(),
    ) else {
        return Ok(None);
    };
    Ok(pdbiox_geom::dihedral(
        previous_carbon,
        nitrogen,
        alpha,
        carbon,
    ))
}

fn psi(
    structure: &Structure,
    current: ResidueRef<'_>,
    next: ResidueRef<'_>,
) -> Result<Option<f64>, RamachandranError> {
    let nitrogen =
        crate::backbone::role_atom(structure, current, PolymerAtomRole::PROTEIN_NITROGEN)
            .map_err(RamachandranError::Roles)?;
    let alpha =
        crate::backbone::role_atom(structure, current, PolymerAtomRole::PROTEIN_ALPHA_CARBON)
            .map_err(RamachandranError::Roles)?;
    let carbon =
        crate::backbone::role_atom(structure, current, PolymerAtomRole::PROTEIN_CARBONYL_CARBON)
            .map_err(RamachandranError::Roles)?;
    let next_nitrogen =
        crate::backbone::role_atom(structure, next, PolymerAtomRole::PROTEIN_NITROGEN)
            .map_err(RamachandranError::Roles)?;
    let (Some(nitrogen), Some(alpha), Some(carbon), Some(next_nitrogen)) =
        (nitrogen, alpha, carbon, next_nitrogen)
    else {
        return Ok(None);
    };
    if !crate::backbone::peptide_bonded(structure, carbon, next_nitrogen) {
        return Ok(None);
    }
    let (Some(nitrogen), Some(alpha), Some(carbon), Some(next_nitrogen)) = (
        nitrogen.position(),
        alpha.position(),
        carbon.position(),
        next_nitrogen.position(),
    ) else {
        return Ok(None);
    };
    Ok(pdbiox_geom::dihedral(
        nitrogen,
        alpha,
        carbon,
        next_nitrogen,
    ))
}

#[cfg(test)]
#[path = "ramachandran_tests.rs"]
mod tests;
