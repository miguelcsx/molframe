//! Canonical voxel coordinates to Cartesian ångströms, sharing cell conventions.

use super::{DensityMap, MrcError, MrcMapDescriptor};
use crate::{CellTransform, numeric::usize_to_f64};
use molframe_core::structure::UnitCell;

impl MrcMapDescriptor {
    /// Returns a row-major homogeneous affine mapping local canonical XYZ voxel
    /// coordinates to Cartesian ångströms. Multiply by a column vector `[x,y,z,1]`.
    ///
    /// File axis permutations have already been canonicalised by the reader.
    /// A nonzero explicit origin replaces the start-index translation, exactly
    /// as in density-map sampling; otherwise starts locate the first voxel.
    ///
    /// # Errors
    ///
    /// Returns an invalid-header error for zero sampling, nonfinite origin, or
    /// a degenerate unit cell, including caller-constructed descriptors.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use molframe_xtal::{MrcBlockReader, MrcBlockOptions};
    /// let reader = MrcBlockReader::open("density.mrc", MrcBlockOptions::default())?;
    /// let affine = reader.descriptor().voxel_to_world()?;
    /// let first_voxel = [affine[0][3], affine[1][3], affine[2][3]];
    /// # Ok::<(), molframe_xtal::MrcError>(())
    /// ```
    pub fn voxel_to_world(&self) -> Result<[[f64; 4]; 4], MrcError> {
        affine(&self.cell, self.sampling, self.starts, self.origin)
    }
}

impl DensityMap {
    /// Returns the row-major local-XYZ voxel to Cartesian ångström affine.
    ///
    /// A nonzero explicit origin takes precedence over start indices. The
    /// matrix uses the same triclinic cell convention as Cartesian sampling.
    ///
    /// # Errors
    ///
    /// Returns an invalid-header error for invalid cell, sampling, or origin.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use molframe_xtal::DensityMap;
    /// let bytes = std::fs::read("density.mrc")?;
    /// let map = DensityMap::from_mrc_bytes(&bytes)?;
    /// let affine = map.voxel_to_world()?;
    /// assert_eq!(affine[3], [0.0, 0.0, 0.0, 1.0]);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn voxel_to_world(&self) -> Result<[[f64; 4]; 4], MrcError> {
        affine(&self.cell, self.sampling, self.starts, self.origin)
    }
}

fn affine(
    cell: &UnitCell,
    sampling: [usize; 3],
    starts: [i32; 3],
    origin: [f64; 3],
) -> Result<[[f64; 4]; 4], MrcError> {
    if sampling.contains(&0) || !origin.iter().all(|value| value.is_finite()) {
        return Err(MrcError::InvalidHeader);
    }
    let transform = CellTransform::new(cell).map_err(|_| MrcError::InvalidHeader)?;
    let inverse_sampling = sampling.map(|count| usize_to_f64(count).recip());
    let translation = if origin.iter().any(|value| value.abs() > f64::EPSILON) {
        origin
    } else {
        transform.to_cartesian(std::array::from_fn(|axis| {
            f64::from(starts[axis]) * inverse_sampling[axis]
        }))
    };
    let mut matrix = [[0.0; 4]; 4];
    for (axis, row) in matrix.iter_mut().take(3).enumerate() {
        for column in 0..3 {
            row[column] = transform.forward_matrix()[axis][column] * inverse_sampling[column];
        }
        row[3] = translation[axis];
    }
    matrix[3][3] = 1.0;
    Ok(matrix)
}

#[cfg(test)]
#[path = "mrc_affine_tests.rs"]
mod tests;
