//! Public TNG trajectory metadata and encoder choices.

use super::TngError;
use crate::Timestep;

/// Compression applied to positions and velocities.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum TngCompression {
    /// Exact IEEE-754 values without container compression.
    Uncompressed,
    /// Exact IEEE-754 values compressed with gzip.
    #[default]
    Lossless,
    /// Native TNG lossy compression with the given inverse-length precision.
    Lossy {
        /// Quantization multiplier in inverse file length units.
        precision: f64,
    },
}

/// Declarative controls for deterministic TNG encoding.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TngWriteOptions {
    /// Base-ten exponent of the file length unit in metres; `-10` is ångström.
    pub distance_unit_exponent: i64,
    /// Position and velocity compression.
    pub compression: TngCompression,
    /// Whether block MD5 hashes are written and verified.
    pub hashes: bool,
}

impl Default for TngWriteOptions {
    fn default() -> Self {
        Self {
            distance_unit_exponent: -9,
            compression: TngCompression::Lossless,
            hashes: true,
        }
    }
}

/// Parsed TNG frames and format metadata not carried by [`Timestep`].
#[derive(Clone, Debug, PartialEq)]
pub struct TngTrajectory {
    /// Coordinate-bearing frames in increasing simulation-step order.
    pub frames: Vec<Timestep>,
    /// Simulation step for each frame.
    pub steps: Vec<i64>,
    /// Base-ten exponent of the source length unit in metres.
    pub distance_unit_exponent: i64,
    /// Lossy coordinate precision recorded in the file header.
    pub compression_precision: f64,
    /// Coordinate block codec retained for faithful rewriting.
    pub compression: TngCompression,
}

/// Multiplier taking a TNG length exponent (metres) to angstrom.
pub(super) fn length_to_angstrom(exponent: i64) -> Result<f64, TngError> {
    let power = i32::try_from(exponent.checked_add(10).ok_or(TngError::InvalidValue)?)
        .map_err(|_| TngError::InvalidValue)?;
    let scale = 10.0_f64.powi(power);
    if scale.is_finite() && scale > 0.0 {
        Ok(scale)
    } else {
        Err(TngError::InvalidValue)
    }
}
