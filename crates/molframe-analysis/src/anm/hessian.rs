//! The matrix-free anisotropic Hessian.
//!
//! Each spring contributes a rank-one `3x3` block along its bond direction, so
//! the whole operator is described by one unit vector per edge. Those vectors
//! are computed once at build time — the square root that normalises them is
//! the expensive part, and an eigensolver applies the operator hundreds of
//! times. Storing them as `f32` halves the bytes the hot loop streams while the
//! accumulation stays in `f64`.
//!
//! One application costs `O(N + E)`: fifteen or so floating-point operations
//! per edge and no indirection beyond the two endpoint lookups.

use crate::network::{ContactGraph, NetworkBudget, NetworkError, check_memory, checked_sum};
use crate::numeric::f64_to_f32;
use crate::{AnmError, network::checked_product};

/// The contact graph plus the bond direction every spring pulls along.
#[derive(Debug)]
pub(super) struct Hessian<'a> {
    pub(super) graph: &'a ContactGraph,
    /// One unit bond vector per edge, in edge order.
    directions: Vec<[f32; 3]>,
}

impl<'a> Hessian<'a> {
    /// Normalises every edge once.
    ///
    /// # Errors
    ///
    /// Returns [`AnmError::CoincidentSites`] when two contacting sites sit at
    /// the same point, which leaves the spring without a direction to pull
    /// along, and a memory error when the direction table will not fit.
    pub(super) fn build(
        graph: &'a ContactGraph,
        positions: &[[f32; 3]],
        selected: &[u32],
        budget: NetworkBudget,
    ) -> Result<Self, AnmError> {
        let direction_bytes = checked_product(&[graph.edges.len(), size_of::<[f32; 3]>()], budget)?;
        check_memory(
            checked_sum(&[graph.owned_bytes(), direction_bytes], budget)?,
            budget,
        )?;

        let mut directions = Vec::new();
        directions
            .try_reserve_exact(graph.edges.len())
            .map_err(|_| NetworkError::MemoryLimit {
                required: direction_bytes,
                limit: budget.memory_limit_bytes,
            })?;
        for edge in &graph.edges {
            let (left, right) = edge.indices();
            let (Some(&first), Some(&second)) = (selected.get(left), selected.get(right)) else {
                return Err(molframe_spatial::SpatialError::NumericRangeExceeded.into());
            };
            let (Some(tail), Some(head)) = (
                site_position(positions, first),
                site_position(positions, second),
            ) else {
                return Err(molframe_spatial::SpatialError::NumericRangeExceeded.into());
            };
            let delta = [
                f64::from(head[0]) - f64::from(tail[0]),
                f64::from(head[1]) - f64::from(tail[1]),
                f64::from(head[2]) - f64::from(tail[2]),
            ];
            let length = delta.iter().map(|value| value * value).sum::<f64>().sqrt();
            if !length.is_finite() || length <= 0.0 {
                return Err(AnmError::CoincidentSites { first, second });
            }
            directions.push([
                f64_to_f32(delta[0] / length),
                f64_to_f32(delta[1] / length),
                f64_to_f32(delta[2] / length),
            ]);
        }
        Ok(Self { graph, directions })
    }

    /// Sites in the network.
    pub(super) fn site_count(&self) -> usize {
        self.graph.site_count()
    }

    /// Degrees of freedom the operator acts on.
    pub(super) fn dimension(&self) -> usize {
        self.site_count() * 3
    }

    /// Bytes the direction table owns beyond the graph.
    pub(super) fn owned_bytes(&self) -> usize {
        self.graph.owned_bytes() + self.directions.capacity() * size_of::<[f32; 3]>()
    }

    /// Accumulates `out += H * rhs` over contiguous `3N` vectors.
    ///
    /// The caller owns `out`'s initial contents, which lets the shifted
    /// operator seed it with `shift * rhs` and subtract this in one pass.
    pub(super) fn accumulate(&self, out: &mut [f64], rhs: &[f64], scale: f64) {
        for (edge, direction) in self.graph.edges.iter().zip(&self.directions) {
            let (left, right) = edge.indices();
            let (tail, head) = (left * 3, right * 3);
            let axis = [
                f64::from(direction[0]),
                f64::from(direction[1]),
                f64::from(direction[2]),
            ];
            let projected = axis[0] * (rhs[tail] - rhs[head])
                + axis[1] * (rhs[tail + 1] - rhs[head + 1])
                + axis[2] * (rhs[tail + 2] - rhs[head + 2]);
            let weighted = projected * scale;
            for (offset, component) in axis.iter().enumerate() {
                let pull = component * weighted;
                out[tail + offset] += pull;
                out[head + offset] -= pull;
            }
        }
    }
}

fn site_position(positions: &[[f32; 3]], site: u32) -> Option<[f32; 3]> {
    usize::try_from(site)
        .ok()
        .and_then(|index| positions.get(index))
        .copied()
}

#[cfg(test)]
#[path = "hessian_tests.rs"]
mod tests;
