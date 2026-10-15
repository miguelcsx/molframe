//! H5MD errors.

/// A malformed, unsupported, or unwritable H5MD trajectory.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum H5mdError {
    /// The HDF5 container could not be decoded or encoded.
    #[error("HDF5 error: {0}")]
    Hdf5(String),
    /// Required H5MD version or creator metadata is invalid.
    #[error("invalid H5MD metadata")]
    InvalidMetadata,
    /// A required dataset is absent.
    #[error("required H5MD dataset '{0}' is absent")]
    MissingDataset(String),
    /// A dataset shape contradicts the H5MD element definition.
    #[error("H5MD dataset '{dataset}' has shape {found:?}, expected {expected}")]
    InvalidShape {
        /// Absolute dataset path.
        dataset: String,
        /// Required shape description.
        expected: &'static str,
        /// Actual dimensions.
        found: Vec<u64>,
    },
    /// Declared units do not match the explicitly requested unit system.
    #[error("H5MD dataset '{dataset}' has unsupported units '{units}'")]
    InvalidUnits {
        /// Absolute dataset path.
        dataset: String,
        /// Declared units, or `<missing>`.
        units: String,
    },
    /// Frame arrays do not share a common atom count or optional-stream schema.
    #[error("H5MD frames do not share one consistent schema")]
    InconsistentFrames,
    /// A group name, unit scale, cell, or number is invalid.
    #[error("H5MD contains an invalid value")]
    InvalidValue,
}

impl From<hdf5_reader::error::Error> for H5mdError {
    fn from(error: hdf5_reader::error::Error) -> Self {
        Self::Hdf5(error.to_string())
    }
}

impl From<hdf5_writer::Error> for H5mdError {
    fn from(error: hdf5_writer::Error) -> Self {
        Self::Hdf5(error.to_string())
    }
}
