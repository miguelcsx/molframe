//! Stable diagnostic codes for every error the interoperability layer returns.

use crate::{ColumnTableError, DatasetError, DlpackError, GraphError};
use molframe_core::{Code, Diagnostic, diagnostic_from};

diagnostic_from!(ColumnTableError, |_error| Code::E5102);
diagnostic_from!(DlpackError, |_error| Code::E5102);
diagnostic_from!(DatasetError, |error| match error {
    DatasetError::Alignment(inner) => Diagnostic::from(inner).code(),
    DatasetError::Io(_)
    | DatasetError::Json(_)
    | DatasetError::InvalidEntry { .. }
    | DatasetError::DuplicateId { .. }
    | DatasetError::MissingSplitMetadata { .. } => Code::E7101,
    DatasetError::IndexOutOfBounds { .. } => Code::E6009,
    DatasetError::InvalidRatios
    | DatasetError::InvalidThreshold
    | DatasetError::InvalidBatchSize
    | DatasetError::InvalidFilter => Code::E5101,
});
diagnostic_from!(GraphError, |error| match error {
    GraphError::InvalidParameter => Code::E5101,
    GraphError::UnsupportedFeature { .. } | GraphError::MissingFeature { .. } => Code::E6103,
    GraphError::MissingCell => Code::E5004,
    GraphError::RaggedCoordinates => Code::E5102,
    GraphError::Spatial(inner) => Diagnostic::from(inner).code(),
    GraphError::IndexOverflow => Code::E1903,
    GraphError::ResourceLimit => Code::E7001,
});

#[cfg(test)]
#[path = "diagnostics_tests.rs"]
mod tests;
