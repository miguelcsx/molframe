//! Exact inverse-image matching for explicit symmetry representatives.

use crate::symmetry::{Rational, SymmetrySet, symmetry_capacity, symmetry_error};
use pdbiox_core::Diagnostic;

impl SymmetrySet {
    pub(crate) fn inverse_image(
        &self,
        operation: usize,
        lattice: [i32; 3],
    ) -> Result<(usize, [i32; 3]), Diagnostic> {
        let source = self
            .operations
            .get(operation)
            .ok_or_else(symmetry_capacity)?;
        let inverse_rotation = inverse_rotation(source.rotation)
            .ok_or_else(|| symmetry_error("symmetry rotation has no integer inverse"))?;
        let mut inverse_translation = [Rational::ZERO; 3];
        for row in 0..3 {
            let mut value = Rational::ZERO;
            for axis in 0..3 {
                let shifted =
                    source.translation[axis].checked_add(Rational::new(lattice[axis], 1)?)?;
                value = value.checked_add(shifted.checked_mul(-inverse_rotation[row][axis])?)?;
            }
            inverse_translation[row] = value;
        }
        for (candidate_index, candidate) in self.operations.iter().enumerate() {
            if candidate.rotation != inverse_rotation {
                continue;
            }
            let mut shift = [0; 3];
            let mut equivalent = true;
            for axis in 0..3 {
                let difference =
                    inverse_translation[axis].checked_sub(candidate.translation[axis])?;
                let Some(integer) = difference.integer() else {
                    equivalent = false;
                    break;
                };
                shift[axis] = integer;
            }
            if equivalent {
                return Ok((candidate_index, shift));
            }
        }
        Err(symmetry_error(
            "explicit symmetry representatives are not closed under inversion",
        ))
    }
}

fn inverse_rotation(matrix: [[i32; 3]; 3]) -> Option<[[i32; 3]; 3]> {
    let determinant = determinant(matrix);
    if determinant.unsigned_abs() != 1 {
        return None;
    }
    let adjugate = [
        [
            matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1],
            matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2],
            matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1],
        ],
        [
            matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2],
            matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0],
            matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2],
        ],
        [
            matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0],
            matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1],
            matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0],
        ],
    ];
    Some(adjugate.map(|row| row.map(|value| value / determinant)))
}

fn determinant(matrix: [[i32; 3]; 3]) -> i32 {
    matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
}
