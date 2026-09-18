//! Dielectric response from dipole fluctuations with an explicit unit prefactor.

use crate::numeric::f64_from_usize;

/// Unit-system and ensemble prefactor for dielectric estimation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DielectricOptions {
    /// Multiplier applied to `���|M|²��� - |���M���|²` before adding one.
    ///
    /// In a chosen unit system this normally contains volume, temperature,
    /// Boltzmann's constant, vacuum permittivity, and the isotropic factor.
    pub fluctuation_prefactor: f64,
}

/// Dipole statistics and isotropic relative permittivity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DielectricEstimate {
    /// Mean total dipole vector.
    pub mean_dipole: [f64; 3],
    /// Population covariance tensor of the dipole components.
    pub covariance: [[f64; 3]; 3],
    /// Scalar fluctuation `trace(covariance)`.
    pub fluctuation: f64,
    /// `1 + fluctuation_prefactor * fluctuation`.
    pub relative_permittivity: f64,
}

/// Why dielectric estimation could not be evaluated.
#[derive(Clone, Copy, Debug, thiserror::Error, PartialEq, Eq)]
pub enum DielectricError {
    /// At least one dipole observation is required.
    #[error("at least one dipole observation is required")]
    Empty,
    /// Dipoles and the prefactor must be finite, and the prefactor non-negative.
    #[error("dipoles and prefactor must be finite and prefactor non-negative")]
    InvalidInput,
}

/// Estimates isotropic relative permittivity from total-dipole observations.
///
/// The physical prefactor is mandatory rather than synthesized from assumed
/// temperature, volume, ensemble, boundary conditions, or units.
///
/// # Errors
///
/// Returns [`DielectricError`] for absent or non-finite data or invalid prefactor.
pub fn dielectric_from_dipoles(
    dipoles: &[[f64; 3]],
    options: DielectricOptions,
) -> Result<DielectricEstimate, DielectricError> {
    if dipoles.is_empty() {
        return Err(DielectricError::Empty);
    }
    if !options.fluctuation_prefactor.is_finite()
        || options.fluctuation_prefactor < 0.0
        || dipoles.iter().flatten().any(|value| !value.is_finite())
    {
        return Err(DielectricError::InvalidInput);
    }
    let count = f64_from_usize(dipoles.len()).ok_or(DielectricError::InvalidInput)?;
    let mean_dipole =
        std::array::from_fn(|axis| dipoles.iter().map(|dipole| dipole[axis]).sum::<f64>() / count);
    let covariance = std::array::from_fn(|row| {
        std::array::from_fn(|column| {
            dipoles
                .iter()
                .map(|dipole| {
                    (dipole[row] - mean_dipole[row]) * (dipole[column] - mean_dipole[column])
                })
                .sum::<f64>()
                / count
        })
    });
    let fluctuation = (0..3).map(|axis| covariance[axis][axis]).sum();
    Ok(DielectricEstimate {
        mean_dipole,
        covariance,
        fluctuation,
        relative_permittivity: 1.0 + options.fluctuation_prefactor * fluctuation,
    })
}

#[cfg(test)]
#[path = "dielectric_tests.rs"]
mod tests;
