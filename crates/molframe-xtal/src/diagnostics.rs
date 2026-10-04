//! Stable diagnostic codes for the map and reflection readers.

use crate::{MapStatisticsError, MrcError, ReflectionError, UnknownColumnType};
use molframe_core::{Code, diagnostic_from};

diagnostic_from!(MrcError, |error| match error {
    MrcError::Truncated | MrcError::Io(_) => Code::E7101,
    MrcError::InvalidHeader => Code::E1201,
    MrcError::UnsupportedMode(_) => Code::E4105,
    MrcError::SizeOverflow => Code::E1903,
    MrcError::ResourceLimit => Code::E1901,
    MrcError::NonFiniteDensity | MrcError::InvalidRegion | MrcError::InvalidMemoryLimit { .. } =>
        Code::E5101,
    MrcError::MemoryLimit { .. } => Code::E7001,
});

diagnostic_from!(MapStatisticsError, |error| match error {
    MapStatisticsError::Empty => Code::E5103,
    MapStatisticsError::MaskLength => Code::E5102,
    MapStatisticsError::NonFiniteDensity | MapStatisticsError::InvalidHistogram => Code::E5101,
});

diagnostic_from!(ReflectionError, |error| match error {
    ReflectionError::MissingTable => Code::E2001,
    ReflectionError::ColumnLength => Code::E5102,
    ReflectionError::MillerIndices => Code::E2004,
    ReflectionError::Metadata => Code::E2002,
    ReflectionError::InvalidNumber => Code::E5101,
    ReflectionError::Unsupported(_) => Code::E4105,
    ReflectionError::InvalidMtz => Code::E1201,
    ReflectionError::TruncatedMtz => Code::E7101,
});

diagnostic_from!(UnknownColumnType, |_error| Code::E5101);

#[cfg(test)]
#[path = "diagnostics_tests.rs"]
mod tests;
