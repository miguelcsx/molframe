//! Residence and survival of explicitly classified water molecules.

use crate::numeric::f64_from_u64;

/// Whether temporary departures break a survival event.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SurvivalMode {
    /// A water must remain present at every frame through the lag.
    Continuous,
    /// A water need only be present at the origin and endpoint.
    Intermittent,
}

/// Survival probability at one lag.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterSurvival {
    /// Frame lag.
    pub lag: usize,
    /// Initially present water-origin observations.
    pub origins: u64,
    /// Probability, or `None` when there are no initially present observations.
    pub probability: Option<f64>,
}

/// Survival series and its trapezoidal time integral.
#[derive(Clone, Debug, PartialEq)]
pub struct WaterDynamics {
    /// Survival values for lags `0..=maximum_lag`.
    pub survival: Vec<WaterSurvival>,
    /// Integrated residence time, absent if any required lag has no observations.
    pub residence_time: Option<f64>,
}

/// Why water dynamics could not be evaluated.
#[derive(Clone, Copy, Debug, thiserror::Error, PartialEq, Eq)]
pub enum WaterDynamicsError {
    /// At least one frame and one explicitly tracked water are required.
    #[error("occupancy must contain at least one frame and one water")]
    Empty,
    /// Every frame must track the same waters.
    #[error("all occupancy frames must have equal length")]
    DimensionMismatch,
    /// Maximum lag must be smaller than the frame count.
    #[error("maximum lag must be smaller than the frame count")]
    InvalidLag,
    /// Frame duration must be finite and positive.
    #[error("frame duration must be finite and positive")]
    InvalidFrameDuration,
    /// Survival counts exceed exact numeric representation.
    #[error("water survival count exceeds exact numeric representation")]
    ObservationLimit,
}

/// Computes water survival and residence time from an explicit occupancy matrix.
///
/// Rows are frames and columns are stable water identities. Water recognition,
/// site selection, and geometric membership are deliberately performed before
/// this kernel, so no solvent names or distance thresholds are embedded here.
///
/// # Errors
///
/// Returns [`WaterDynamicsError`] for empty/ragged input, invalid lag, or time step.
pub fn water_dynamics(
    occupancy: &[Vec<bool>],
    maximum_lag: usize,
    frame_duration: f64,
    mode: SurvivalMode,
) -> Result<WaterDynamics, WaterDynamicsError> {
    validate(occupancy, maximum_lag, frame_duration)?;
    let survival: Vec<WaterSurvival> = (0..=maximum_lag)
        .map(|lag| survival_at_lag(occupancy, lag, mode))
        .collect::<Result<_, _>>()?;
    let residence_time = integrate(&survival, frame_duration);
    Ok(WaterDynamics {
        survival,
        residence_time,
    })
}

fn validate(
    occupancy: &[Vec<bool>],
    maximum_lag: usize,
    frame_duration: f64,
) -> Result<(), WaterDynamicsError> {
    let Some(first) = occupancy.first() else {
        return Err(WaterDynamicsError::Empty);
    };
    if first.is_empty() {
        return Err(WaterDynamicsError::Empty);
    }
    if occupancy.iter().any(|frame| frame.len() != first.len()) {
        return Err(WaterDynamicsError::DimensionMismatch);
    }
    if maximum_lag >= occupancy.len() {
        return Err(WaterDynamicsError::InvalidLag);
    }
    if !frame_duration.is_finite() || frame_duration <= 0.0 {
        return Err(WaterDynamicsError::InvalidFrameDuration);
    }
    Ok(())
}

fn survival_at_lag(
    occupancy: &[Vec<bool>],
    lag: usize,
    mode: SurvivalMode,
) -> Result<WaterSurvival, WaterDynamicsError> {
    let mut origins = 0_u64;
    let mut surviving = 0_u64;
    for origin in 0..occupancy.len() - lag {
        for water in 0..occupancy[origin].len() {
            if !occupancy[origin][water] {
                continue;
            }
            origins = origins
                .checked_add(1)
                .ok_or(WaterDynamicsError::ObservationLimit)?;
            let survives = match mode {
                SurvivalMode::Continuous => occupancy[origin..=origin + lag]
                    .iter()
                    .all(|frame| frame[water]),
                SurvivalMode::Intermittent => occupancy[origin + lag][water],
            };
            surviving = surviving
                .checked_add(u64::from(survives))
                .ok_or(WaterDynamicsError::ObservationLimit)?;
        }
    }
    let probability = if origins == 0 {
        None
    } else {
        let surviving = f64_from_u64(surviving).ok_or(WaterDynamicsError::ObservationLimit)?;
        let origins = f64_from_u64(origins).ok_or(WaterDynamicsError::ObservationLimit)?;
        Some(surviving / origins)
    };
    Ok(WaterSurvival {
        lag,
        origins,
        probability,
    })
}

fn integrate(survival: &[WaterSurvival], frame_duration: f64) -> Option<f64> {
    survival.windows(2).try_fold(0.0, |area, pair| {
        let left = pair[0].probability?;
        let right = pair[1].probability?;
        Some(area + frame_duration * (left + right) / 2.0)
    })
}

#[cfg(test)]
#[path = "water_dynamics_tests.rs"]
mod tests;
