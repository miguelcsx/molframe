//! Direct X-ray structure-factor summation over a model and its space group.

use crate::{CellTransform, GaussianFormFactor, SymmetryOperation};
use molframe_core::Element;
use std::f64::consts::PI;

/// A complex amplitude.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Complex64 {
    /// Real part.
    pub re: f64,
    /// Imaginary part.
    pub im: f64,
}

impl Complex64 {
    /// The modulus `|F|`.
    #[must_use]
    pub fn norm(self) -> f64 {
        self.re.hypot(self.im)
    }

    /// The phase angle in radians, in `(−π, π]`.
    #[must_use]
    pub fn arg(self) -> f64 {
        self.im.atan2(self.re)
    }
}

/// How an atom's thermal motion damps its scattering.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Displacement {
    /// Isotropic `B` in square ångström.
    Isotropic(f64),
    /// Anisotropic `U` in orthogonal square ångström, ordered
    /// `[U11, U22, U33, U12, U13, U23]`.
    Anisotropic([f64; 6]),
}

/// One atom of the asymmetric unit as the summation sees it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScatteringSite {
    /// The element, which selects the form factor.
    pub element: Element,
    /// Fractional coordinates.
    pub position: [f64; 3],
    /// Site occupancy.
    pub occupancy: f64,
    /// Thermal displacement.
    pub displacement: Displacement,
}

/// Why a structure factor could not be computed.
#[derive(Clone, Copy, Debug, thiserror::Error, PartialEq, Eq)]
pub enum StructureFactorError {
    /// An element has no tabulated form factor.
    #[error("no X-ray form factor is tabulated for atomic number {0}")]
    MissingFormFactor(u8),
    /// No symmetry operation was supplied.
    #[error("structure factors need the complete set of symmetry operations")]
    NoOperations,
}

/// Sums atomic scattering over a cell's symmetry-equivalent positions.
///
/// `F(hkl) = Σ_sites Σ_ops occ · f(s) · T · exp(2πi h·(W x + w))`, where the
/// operations are the complete space group including the identity and the
/// centring translations, `s = sin θ / λ` and `T` is the Debye–Waller factor of
/// the displacement carried through each operation.
#[derive(Clone, Debug)]
pub struct StructureFactorCalculator<'a> {
    cell: CellTransform,
    operations: &'a [SymmetryOperation],
}

impl<'a> StructureFactorCalculator<'a> {
    /// Binds a cell to the complete operation set of its space group.
    ///
    /// # Errors
    ///
    /// Returns [`StructureFactorError::NoOperations`] for an empty set.
    pub fn new(
        cell: CellTransform,
        operations: &'a [SymmetryOperation],
    ) -> Result<Self, StructureFactorError> {
        if operations.is_empty() {
            return Err(StructureFactorError::NoOperations);
        }
        Ok(Self { cell, operations })
    }

    /// The structure factor of one reflection.
    ///
    /// # Errors
    ///
    /// Returns an error for an element without a tabulated form factor.
    pub fn calculate(
        &self,
        sites: &[ScatteringSite],
        hkl: [i32; 3],
    ) -> Result<Complex64, StructureFactorError> {
        let stol2 = 0.25 * self.cell.reciprocal_spacing_squared(hkl);
        let indices = hkl.map(f64::from);
        let mut total = Complex64::default();
        for site in sites {
            let form = GaussianFormFactor::xray(site.element).ok_or(
                StructureFactorError::MissingFormFactor(site.element.atomic_number()),
            )?;
            let weight = site.occupancy * form.value(stol2);
            let (mut real, mut imaginary) = (0.0, 0.0);
            for operation in self.operations {
                let position = operation.apply_fractional(site.position);
                let angle = 2.0
                    * PI
                    * (indices[0] * position[0]
                        + indices[1] * position[1]
                        + indices[2] * position[2]);
                let damping = self.damping(site.displacement, stol2, operation, indices);
                real += damping * angle.cos();
                imaginary += damping * angle.sin();
            }
            total.re += weight * real;
            total.im += weight * imaginary;
        }
        Ok(total)
    }

    /// The structure factors of many reflections, in input order.
    ///
    /// # Errors
    ///
    /// Returns the first error any reflection raises.
    pub fn calculate_many(
        &self,
        sites: &[ScatteringSite],
        reflections: &[[i32; 3]],
    ) -> Result<Vec<Complex64>, StructureFactorError> {
        reflections
            .iter()
            .map(|hkl| self.calculate(sites, *hkl))
            .collect()
    }

    /// Debye–Waller factor of one site seen through one operation.
    fn damping(
        &self,
        displacement: Displacement,
        stol2: f64,
        operation: &SymmetryOperation,
        indices: [f64; 3],
    ) -> f64 {
        match displacement {
            Displacement::Isotropic(b) => (-b * stol2).exp(),
            Displacement::Anisotropic(u) => {
                // The index vector the operation sends `h` to is `Wᵀ h`; the
                // tensor is carried into fractional axes by `F U Fᵀ`.
                let rotated: [f64; 3] = std::array::from_fn(|column| {
                    (0..3)
                        .map(|row| f64::from(operation.rotation[row][column]) * indices[row])
                        .sum()
                });
                let f = self.cell.inverse_matrix();
                let tensor = [[u[0], u[3], u[4]], [u[3], u[1], u[5]], [u[4], u[5], u[2]]];
                // Fᵀ h' first, then the quadratic form in U.
                let projected: [f64; 3] =
                    std::array::from_fn(|axis| (0..3).map(|row| f[row][axis] * rotated[row]).sum());
                let quadratic: f64 = (0..3)
                    .map(|i| {
                        (0..3)
                            .map(|j| projected[i] * tensor[i][j] * projected[j])
                            .sum::<f64>()
                    })
                    .sum();
                (-2.0 * PI * PI * quadratic).exp()
            }
        }
    }
}

#[cfg(test)]
#[path = "structure_factor_tests.rs"]
mod tests;
