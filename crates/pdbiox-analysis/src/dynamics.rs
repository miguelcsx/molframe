//! Water occupancy dynamics and dipole-fluctuation dielectric analysis.

use pdbiox_core::AtomSelection;

use crate::numeric::{u64_to_f64, usize_to_f64};

/// One lag of stable-ID water occupancy dynamics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterLag {
    /// Frame lag.
    pub lag: usize,
    /// Origin-water observations in the denominator.
    pub observations: u64,
    /// Waters present at both endpoints, allowing intermittent absence.
    pub surviving: u64,
    /// Waters continuously present through the full interval.
    pub resident: u64,
    /// `surviving / observations`, absent for an empty denominator.
    pub survival_probability: Option<f64>,
    /// `resident / observations`, absent for an empty denominator.
    pub residence_probability: Option<f64>,
}

/// Explicit lag range for water dynamics.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WaterDynamicsOptions {
    /// Largest frame lag, including zero.
    pub maximum_lag: usize,
}

/// Dipole fluctuation and resulting relative dielectric constant.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DielectricResult {
    /// Mean dipole vector in caller units.
    pub mean_dipole: [f64; 3],
    /// `<|M|²> - |<M>|²` in squared caller dipole units.
    pub fluctuation: f64,
    /// Relative dielectric constant under the explicit conversion prefactor.
    pub relative_permittivity: f64,
}

/// Explicit unit and thermodynamic definition for dielectric estimation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DielectricOptions {
    /// Mean sampled volume in caller length units cubed.
    pub volume: f64,
    /// Absolute temperature.
    pub temperature: f64,
    /// Unit-system and boundary conversion multiplying fluctuation/(V T).
    pub fluctuation_prefactor: f64,
}

/// Invalid explicit dynamics input.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum DynamicsError {
    /// Frame series or lag range is empty/incompatible.
    #[error("water occupancy frames must be non-empty and cover the requested lag")]
    InvalidLag,
    /// Dipoles or thermodynamic/unit controls are invalid.
    #[error("dielectric inputs must be finite with positive volume, temperature and prefactor")]
    InvalidDielectricInput,
}

/// Computes intermittent survival and continuous residence probabilities.
///
/// Each selection contains stable caller-defined water IDs inside the region at
/// one frame. Chemical water recognition therefore stays outside this kernel.
///
/// # Errors
///
/// Returns [`DynamicsError::InvalidLag`] for an empty series or excessive lag.
pub fn water_dynamics(
    occupancy: &[AtomSelection],
    options: WaterDynamicsOptions,
) -> Result<Vec<WaterLag>, DynamicsError> {
    if occupancy.is_empty() || options.maximum_lag >= occupancy.len() {
        return Err(DynamicsError::InvalidLag);
    }
    Ok((0..=options.maximum_lag)
        .map(|lag| water_lag(occupancy, lag))
        .collect())
}

fn water_lag(occupancy: &[AtomSelection], lag: usize) -> WaterLag {
    let mut observations = 0u64;
    let mut surviving = 0u64;
    let mut resident = 0u64;
    for origin in 0..occupancy.len() - lag {
        for water in &occupancy[origin] {
            observations += 1;
            if occupancy[origin + lag].contains(water) {
                surviving += 1;
            }
            if occupancy[origin..=origin + lag]
                .iter()
                .all(|frame| frame.contains(water))
            {
                resident += 1;
            }
        }
    }
    WaterLag {
        lag,
        observations,
        surviving,
        resident,
        survival_probability: ratio(surviving, observations),
        residence_probability: ratio(resident, observations),
    }
}

fn ratio(numerator: u64, denominator: u64) -> Option<f64> {
    (denominator > 0).then(|| u64_to_f64(numerator) / u64_to_f64(denominator))
}

/// Computes dielectric response from an explicit dipole time series.
///
/// # Errors
///
/// Returns invalid-input errors for empty/non-finite data or non-positive
/// thermodynamic and unit-conversion controls.
pub fn dielectric_from_dipoles(
    dipoles: &[[f64; 3]],
    options: DielectricOptions,
) -> Result<DielectricResult, DynamicsError> {
    if dipoles.is_empty()
        || dipoles.iter().flatten().any(|value| !value.is_finite())
        || !options.volume.is_finite()
        || options.volume <= 0.0
        || !options.temperature.is_finite()
        || options.temperature <= 0.0
        || !options.fluctuation_prefactor.is_finite()
        || options.fluctuation_prefactor <= 0.0
    {
        return Err(DynamicsError::InvalidDielectricInput);
    }
    let count = usize_to_f64(dipoles.len());
    let mean_dipole =
        std::array::from_fn(|axis| dipoles.iter().map(|dipole| dipole[axis]).sum::<f64>() / count);
    let mean_squared = dipoles
        .iter()
        .map(|dipole| dipole.iter().map(|value| value * value).sum::<f64>())
        .sum::<f64>()
        / count;
    let squared_mean = mean_dipole.iter().map(|value| value * value).sum::<f64>();
    let fluctuation = (mean_squared - squared_mean).max(0.0);
    let relative_permittivity =
        1.0 + options.fluctuation_prefactor * fluctuation / (options.volume * options.temperature);
    Ok(DielectricResult {
        mean_dipole,
        fluctuation,
        relative_permittivity,
    })
}

#[cfg(test)]
#[path = "dynamics_tests.rs"]
mod tests;
