//! Public, validated spatial planning profiles.

use crate::{SpatialBackend, SpatialError, SpatialOption};

/// Tunable decisions used when [`SpatialBackend::Auto`] is requested.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct AutoBackendProfile {
    /// Largest Cartesian product evaluated directly.
    pub brute_force_pair_limit: usize,
    /// Smallest target set for which a k-d tree may be selected.
    pub kd_target_minimum: usize,
    /// Required target-to-query size ratio for selecting a k-d tree.
    pub kd_query_ratio: usize,
    /// Backend selected for periodic searches.
    pub periodic_backend: SpatialBackend,
}

impl AutoBackendProfile {
    /// General-purpose profile used by [`Default`].
    pub const BALANCED: Self = Self {
        brute_force_pair_limit: 250_000,
        kd_target_minimum: 20_000,
        kd_query_ratio: 8,
        periodic_backend: SpatialBackend::BruteForce,
    };

    /// Validates that every threshold defines an unambiguous plan.
    ///
    /// # Errors
    ///
    /// Returns the first invalid automatic-planning field.
    pub fn validate(self) -> Result<(), SpatialError> {
        if self.brute_force_pair_limit == 0 {
            return Err(SpatialError::InvalidOption(
                SpatialOption::BruteForcePairLimit,
            ));
        }
        if self.kd_target_minimum == 0 {
            return Err(SpatialError::InvalidOption(SpatialOption::KdTargetMinimum));
        }
        if self.kd_query_ratio == 0 {
            return Err(SpatialError::InvalidOption(SpatialOption::KdQueryRatio));
        }
        if self.periodic_backend == SpatialBackend::Auto {
            return Err(SpatialError::InvalidOption(SpatialOption::PeriodicBackend));
        }
        Ok(())
    }

    /// Resolves one workload after validating the profile.
    ///
    /// # Errors
    ///
    /// Returns an error when the profile is invalid.
    pub fn resolve(
        self,
        left_count: usize,
        right_count: usize,
        periodic: bool,
    ) -> Result<SpatialBackend, SpatialError> {
        self.validate()?;
        if periodic {
            return Ok(self.periodic_backend);
        }
        let pair_count = left_count
            .checked_mul(right_count)
            .ok_or(SpatialError::NumericRangeExceeded)?;
        if pair_count <= self.brute_force_pair_limit {
            return Ok(SpatialBackend::BruteForce);
        }
        if right_count >= self.kd_target_minimum && left_count < right_count / self.kd_query_ratio {
            return Ok(SpatialBackend::KdTree);
        }
        Ok(SpatialBackend::CellList)
    }
}

impl Default for AutoBackendProfile {
    fn default() -> Self {
        Self::BALANCED
    }
}

/// Rule for deriving a one-shot neighbour-list skin from the query cutoff.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct NeighborSkinProfile {
    /// Fraction of the cutoff used as the preferred skin.
    pub cutoff_ratio: f32,
    /// Lower bound for the derived skin, in ångströms.
    pub minimum: f32,
}

impl NeighborSkinProfile {
    /// General-purpose profile used by [`Default`].
    pub const BALANCED: Self = Self {
        cutoff_ratio: 0.2,
        minimum: 0.5,
    };

    /// Validates finite, non-negative skin parameters.
    ///
    /// # Errors
    ///
    /// Returns the first invalid skin field.
    pub fn validate(self) -> Result<(), SpatialError> {
        if !self.cutoff_ratio.is_finite() || self.cutoff_ratio < 0.0 {
            return Err(SpatialError::InvalidOption(
                SpatialOption::NeighborSkinRatio,
            ));
        }
        if !self.minimum.is_finite() || self.minimum < 0.0 {
            return Err(SpatialError::InvalidOption(
                SpatialOption::NeighborSkinMinimum,
            ));
        }
        Ok(())
    }

    /// Computes the skin associated with `cutoff`.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid profile fields or a non-finite result.
    pub fn skin(self, cutoff: f32) -> Result<f32, SpatialError> {
        self.validate()?;
        let skin = (cutoff * self.cutoff_ratio).max(self.minimum);
        if skin.is_finite() {
            Ok(skin)
        } else {
            Err(SpatialError::InvalidOption(
                SpatialOption::NeighborSkinRatio,
            ))
        }
    }
}

impl Default for NeighborSkinProfile {
    fn default() -> Self {
        Self::BALANCED
    }
}

/// Memory and coarsening policy for contiguous cell grids.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct CellGridOptions {
    /// Maximum number of allocated grid cells.
    pub maximum_cell_count: usize,
    /// Edge multiplier used when a grid exceeds the cell budget.
    pub edge_growth_factor: f32,
}

impl CellGridOptions {
    /// General-purpose profile used by [`Default`].
    pub const MEMORY_BALANCED: Self = Self {
        maximum_cell_count: 1_000_000,
        edge_growth_factor: 2.0,
    };

    /// Validates a non-empty budget and a convergent coarsening factor.
    ///
    /// # Errors
    ///
    /// Returns the first invalid grid field.
    pub fn validate(self) -> Result<(), SpatialError> {
        if self.maximum_cell_count == 0 {
            return Err(SpatialError::InvalidOption(SpatialOption::MaximumCellCount));
        }
        if !self.edge_growth_factor.is_finite() || self.edge_growth_factor <= 1.0 {
            return Err(SpatialError::InvalidOption(SpatialOption::CellGrowthFactor));
        }
        Ok(())
    }
}

impl Default for CellGridOptions {
    fn default() -> Self {
        Self::MEMORY_BALANCED
    }
}

/// Image-search budget for periodic k-d tree queries.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct KdPeriodicOptions {
    /// Maximum query images evaluated for one periodic tree query.
    pub maximum_image_count: usize,
}

impl KdPeriodicOptions {
    /// General-purpose bounded image policy used by [`Default`].
    pub const BALANCED: Self = Self {
        maximum_image_count: 4_096,
    };

    /// Validates a non-empty image budget.
    ///
    /// # Errors
    ///
    /// Returns an error when no image can be evaluated.
    pub fn validate(self) -> Result<(), SpatialError> {
        if self.maximum_image_count == 0 {
            Err(SpatialError::InvalidOption(
                SpatialOption::KdPeriodicImageLimit,
            ))
        } else {
            Ok(())
        }
    }
}

impl Default for KdPeriodicOptions {
    fn default() -> Self {
        Self::BALANCED
    }
}

/// Complete construction options for a reusable neighbour list.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct NeighborListOptions {
    /// Absolute displacement skin, in ångströms.
    pub skin: f32,
    /// Cell-grid policy used to construct the candidate set.
    pub cell_grid: CellGridOptions,
}

impl NeighborListOptions {
    /// Creates options with an explicit skin and the default cell-grid policy.
    #[must_use]
    pub const fn with_skin(skin: f32) -> Self {
        Self {
            skin,
            cell_grid: CellGridOptions::MEMORY_BALANCED,
        }
    }

    /// Validates the skin and cell-grid fields.
    ///
    /// # Errors
    ///
    /// Returns the first invalid field.
    pub fn validate(self) -> Result<(), SpatialError> {
        if !self.skin.is_finite() || self.skin < 0.0 {
            return Err(SpatialError::InvalidOption(
                SpatialOption::NeighborSkinMinimum,
            ));
        }
        self.cell_grid.validate()
    }
}

/// Complete options for a one-shot spatial query.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct SpatialSearchOptions {
    /// Requested implementation, including inspectable automatic selection.
    pub backend: SpatialBackend,
    /// Automatic backend thresholds.
    pub automatic: AutoBackendProfile,
    /// One-shot neighbour-list skin rule.
    pub neighbor_skin: NeighborSkinProfile,
    /// Cell-grid memory and coarsening policy.
    pub cell_grid: CellGridOptions,
    /// Periodic k-d tree image-search budget.
    pub kd_periodic: KdPeriodicOptions,
}

impl SpatialSearchOptions {
    /// General-purpose automatic search profile used by [`Default`].
    pub const BALANCED: Self = Self {
        backend: SpatialBackend::Auto,
        automatic: AutoBackendProfile::BALANCED,
        neighbor_skin: NeighborSkinProfile::BALANCED,
        cell_grid: CellGridOptions::MEMORY_BALANCED,
        kd_periodic: KdPeriodicOptions::BALANCED,
    };

    /// Returns the balanced profile with an explicit backend.
    #[must_use]
    pub const fn with_backend(backend: SpatialBackend) -> Self {
        Self {
            backend,
            ..Self::BALANCED
        }
    }

    /// Validates all query policy, including branches unused by this workload.
    ///
    /// # Errors
    ///
    /// Returns the first invalid field.
    pub fn validate(self) -> Result<(), SpatialError> {
        self.automatic.validate()?;
        self.neighbor_skin.validate()?;
        self.cell_grid.validate()?;
        self.kd_periodic.validate()
    }

    /// Produces an inspectable execution plan for one workload shape.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid profile fields or derived values.
    pub fn plan(
        self,
        left_count: usize,
        right_count: usize,
        periodic: bool,
        cutoff: f32,
    ) -> Result<SpatialPlan, SpatialError> {
        if !cutoff.is_finite() || cutoff < 0.0 {
            return Err(SpatialError::InvalidCutoff);
        }
        self.validate()?;
        let backend = match self.backend {
            SpatialBackend::Auto => self.automatic.resolve(left_count, right_count, periodic)?,
            backend => backend,
        };
        let neighbor_skin = if backend == SpatialBackend::NeighborList {
            Some(self.neighbor_skin.skin(cutoff)?)
        } else {
            None
        };
        Ok(SpatialPlan {
            requested_backend: self.backend,
            backend,
            neighbor_skin,
            options: self,
        })
    }
}

impl Default for SpatialSearchOptions {
    fn default() -> Self {
        Self::BALANCED
    }
}

/// Resolved, provenance-friendly decision for a spatial workload.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct SpatialPlan {
    /// Backend requested by the caller.
    pub requested_backend: SpatialBackend,
    /// Concrete backend selected for execution.
    pub backend: SpatialBackend,
    /// Derived neighbour skin when the selected backend requires one.
    pub neighbor_skin: Option<f32>,
    /// Complete validated policy that produced this decision.
    pub options: SpatialSearchOptions,
}

#[cfg(test)]
#[path = "options_tests.rs"]
mod tests;
