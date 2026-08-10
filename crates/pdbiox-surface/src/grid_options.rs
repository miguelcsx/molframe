//! Explicit resolution and allocation limits for voxel surface algorithms.

/// Numerical resolution and hard allocation ceiling for a surface grid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceGridOptions {
    /// Cubic grid-cell edge length in ångström.
    pub resolution: f32,
    /// Maximum number of cells the operation may allocate.
    pub max_cells: usize,
}

impl SurfaceGridOptions {
    /// Cell ceiling used by the named [`Self::standard`] profile.
    pub const STANDARD_MAX_CELLS: usize = 8_000_000;

    /// Conventional bounded profile used by the convenience surface wrappers.
    #[must_use]
    pub const fn standard(resolution: f32) -> Self {
        Self {
            resolution,
            max_cells: Self::STANDARD_MAX_CELLS,
        }
    }

    pub(crate) fn validate(self) -> Result<Self, crate::SasaError> {
        if self.resolution.is_finite() && self.resolution > 0.0 && self.max_cells > 0 {
            Ok(self)
        } else {
            Err(crate::SasaError::InvalidGridOptions)
        }
    }
}
