//! Stable diagnostic codes for the geometry errors that had none.

use crate::{BatchGeometryError, FluctuationError, PeriodicError, RotationError};
use molframe_core::{Code, diagnostic_from};
use std::fmt;

impl fmt::Display for FluctuationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::NoFrames => "no frames were supplied, so there is nothing to average over",
            Self::RaggedFrames => "frames disagree on how many atoms they hold",
            Self::TooManyFrames => "the frame count exceeds what the accumulator can represent",
        })
    }
}

impl fmt::Display for RotationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::NonFinite => "a rotation contained NaN or infinity",
            Self::NotOrthonormal => "the matrix is not orthonormal within the tolerance",
            Self::Reflection => "the matrix is a reflection, not a rotation",
            Self::DidNotConverge => "the rotation mean did not converge",
            Self::Empty => "no rotations were supplied",
            Self::InvalidOptions => "the rotation mean options are invalid",
            Self::TooManyRotations => "the rotation count exceeds the numeric kernel's range",
        })
    }
}

impl fmt::Display for PeriodicError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::NonFinite => "an angle, tolerance or weight was not finite",
            Self::DimensionMismatch => "points or a metric have different dimensions",
            Self::InvalidWeight => "a metric weight was zero or negative",
            Self::Empty => "no observation was supplied",
            Self::TooManyObservations => "the observation count exceeds the numeric kernel's range",
        })
    }
}

diagnostic_from!(BatchGeometryError, |_error| Code::E5102);
diagnostic_from!(FluctuationError, |error| match error {
    FluctuationError::NoFrames => Code::E5103,
    FluctuationError::RaggedFrames => Code::E5102,
    FluctuationError::TooManyFrames => Code::E1903,
});
diagnostic_from!(RotationError, |error| match error {
    RotationError::NonFinite
    | RotationError::NotOrthonormal
    | RotationError::Reflection
    | RotationError::InvalidOptions => Code::E5101,
    RotationError::TooManyRotations => Code::E1903,
    RotationError::DidNotConverge => Code::E5104,
    RotationError::Empty => Code::E5103,
});
diagnostic_from!(PeriodicError, |error| match error {
    PeriodicError::NonFinite | PeriodicError::InvalidWeight => Code::E5101,
    PeriodicError::DimensionMismatch => Code::E5102,
    PeriodicError::Empty => Code::E5103,
    PeriodicError::TooManyObservations => Code::E1903,
});

#[cfg(test)]
#[path = "diagnostics_tests.rs"]
mod tests;
