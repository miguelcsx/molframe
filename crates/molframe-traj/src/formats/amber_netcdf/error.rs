//! AMBER `NetCDF` errors.

/// A malformed, unsupported, or unwritable AMBER `NetCDF` trajectory.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum AmberNetcdfError {
    /// The `NetCDF` container could not be decoded or encoded.
    #[error("NetCDF error: {0}")]
    Netcdf(String),
    /// Required AMBER convention metadata is absent or invalid.
    #[error("invalid AMBER NetCDF convention metadata")]
    InvalidConvention,
    /// A required variable is absent.
    #[error("AMBER NetCDF variable '{0}' is required")]
    MissingVariable(&'static str),
    /// A variable has dimensions that contradict the AMBER convention.
    #[error("AMBER NetCDF variable '{variable}' has shape {found:?}, expected {expected}")]
    InvalidShape {
        /// Variable name.
        variable: &'static str,
        /// Required dimension description.
        expected: &'static str,
        /// Actual dimensions.
        found: Vec<u64>,
    },
    /// A variable carries an unsupported or absent unit declaration.
    #[error("AMBER NetCDF variable '{variable}' has unsupported units '{units}'")]
    InvalidUnits {
        /// Variable name.
        variable: &'static str,
        /// Declared units, or `<missing>`.
        units: String,
    },
    /// Atom counts or optional arrays differ between frames.
    #[error("AMBER NetCDF frames do not share one consistent schema")]
    InconsistentFrames,
    /// A numeric value is non-finite or cannot be represented safely.
    #[error("AMBER NetCDF contains an invalid numeric value")]
    InvalidValue,
}

impl From<netcdf_reader::Error> for AmberNetcdfError {
    fn from(error: netcdf_reader::Error) -> Self {
        Self::Netcdf(error.to_string())
    }
}

impl From<netcdf_writer::Error> for AmberNetcdfError {
    fn from(error: netcdf_writer::Error) -> Self {
        Self::Netcdf(error.to_string())
    }
}
