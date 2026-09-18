//! Explicit resolution and allocation limits for voxel surface algorithms.

/// Numerical resolution and hard cell ceiling for a surface grid.
///
/// Voxel algorithms additionally enforce a complete 500 MB peak-byte ceiling
/// over resident inputs, temporary buffers and returned geometry. Increasing
/// `max_cells` can never bypass that operation-memory guard.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceGridOptions {
    /// Cubic grid-cell edge length in ångström.
    pub resolution: f32,
    /// Maximum number of cells the operation may allocate.
    pub max_cells: usize,
    /// Maximum accounted peak bytes the operation may allocate.
    ///
    /// A caller who has provisioned the machine for a large map raises this;
    /// no constant in this crate can know what their working set may be.
    pub max_workspace_bytes: usize,
}

impl SurfaceGridOptions {
    /// Cell ceiling used by the named [`Self::standard`] profile.
    pub const STANDARD_MAX_CELLS: usize = 8_000_000;
    /// Peak-byte allowance used by the named [`Self::standard`] profile.
    pub const STANDARD_WORKSPACE_BYTES: usize = 500_000_000;

    /// Conventional bounded profile used by the convenience surface wrappers.
    #[must_use]
    pub const fn standard(resolution: f32) -> Self {
        Self {
            resolution,
            max_cells: Self::STANDARD_MAX_CELLS,
            max_workspace_bytes: Self::STANDARD_WORKSPACE_BYTES,
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
