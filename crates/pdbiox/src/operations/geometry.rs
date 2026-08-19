//! Typed coordinate and ensemble geometry requests for the facade plan.

use super::requests::{CoordinateInput, ExecutionPlanError, PlanInput, ScalarInput};

/// A reusable geometry request over borrowed plan inputs.
#[derive(Clone, Debug)]
pub enum GeometryRequest {
    /// Unweighted geometric centre.
    Centroid {
        /// Coordinate-array slot.
        positions: usize,
    },
    /// Geometric or mass-weighted centre.
    CentreOfMass {
        /// Coordinate-array slot.
        positions: usize,
        /// Optional scalar-array slot.
        masses: Option<usize>,
    },
    /// Geometric or mass-weighted radius of gyration.
    RadiusOfGyration {
        /// Coordinate-array slot.
        positions: usize,
        /// Optional scalar-array slot.
        masses: Option<usize>,
    },
    /// Geometric or mass-weighted inertia tensor.
    InertiaTensor {
        /// Coordinate-array slot.
        positions: usize,
        /// Optional scalar-array slot.
        masses: Option<usize>,
    },
    /// Principal axes of the inertia tensor.
    PrincipalAxes {
        /// Coordinate-array slot.
        positions: usize,
        /// Optional scalar-array slot.
        masses: Option<usize>,
        /// Bounded eigensolver controls.
        options: pdbiox_geom::EigenOptions,
    },
    /// Shape asphericity from the gyration tensor.
    Asphericity {
        /// Coordinate-array slot.
        positions: usize,
        /// Bounded eigensolver controls.
        options: pdbiox_geom::EigenOptions,
    },
    /// Eigenvectors and eigenvalues of the gyration tensor.
    GyrationAxes {
        /// Coordinate-array slot.
        positions: usize,
        /// Bounded eigensolver controls.
        options: pdbiox_geom::EigenOptions,
    },
    /// Dense pairwise distances within one coordinate set.
    DistanceMatrix {
        /// Coordinate-array slot.
        positions: usize,
    },
    /// Dense pairwise distances between two coordinate sets.
    DistanceMatrixBetween {
        /// Left coordinate-array slot.
        left: usize,
        /// Right coordinate-array slot.
        right: usize,
    },
    /// Per-atom RMS fluctuation over corresponding frames.
    Rmsf {
        /// Frame-array slot.
        frames: usize,
    },
}

/// Typed result values emitted by [`GeometryRequest`].
#[derive(Clone, Debug)]
pub enum GeometryValue {
    /// Geometric centre.
    Centroid(Option<[f64; 3]>),
    /// Geometric or mass-weighted centre.
    CentreOfMass(Option<[f64; 3]>),
    /// Radius of gyration.
    RadiusOfGyration(Option<f64>),
    /// Inertia tensor.
    InertiaTensor(Option<[[f64; 3]; 3]>),
    /// Principal axes.
    PrincipalAxes(Option<pdbiox_geom::Decomposition<3>>),
    /// Asphericity.
    Asphericity(Option<f64>),
    /// Gyration axes.
    GyrationAxes(Option<pdbiox_geom::Decomposition<3>>),
    /// Dense distance matrix.
    DistanceMatrix(pdbiox_geom::DistanceMatrix),
    /// Per-atom RMS fluctuation.
    Rmsf(Vec<f64>),
}

pub(super) fn execute(
    operation: &str,
    request: &GeometryRequest,
    input: PlanInput<'_>,
) -> Result<GeometryValue, ExecutionPlanError> {
    match request {
        GeometryRequest::Centroid { positions } => Ok(GeometryValue::Centroid(
            pdbiox_geom::centroid(coordinates(operation, input.arrays, *positions)?),
        )),
        GeometryRequest::CentreOfMass { positions, masses } => {
            let positions = coordinates(operation, input.arrays, *positions)?;
            let masses = scalar_values(operation, input.scalars, *masses, positions.len())?;
            Ok(GeometryValue::CentreOfMass(pdbiox_geom::centre_of_mass(
                positions, masses,
            )))
        }
        GeometryRequest::RadiusOfGyration { positions, masses } => {
            let positions = coordinates(operation, input.arrays, *positions)?;
            let masses = scalar_values(operation, input.scalars, *masses, positions.len())?;
            Ok(GeometryValue::RadiusOfGyration(
                pdbiox_geom::radius_of_gyration(positions, masses),
            ))
        }
        GeometryRequest::InertiaTensor { positions, masses } => {
            let positions = coordinates(operation, input.arrays, *positions)?;
            let masses = scalar_values(operation, input.scalars, *masses, positions.len())?;
            Ok(GeometryValue::InertiaTensor(pdbiox_geom::inertia_tensor(
                positions, masses,
            )))
        }
        GeometryRequest::PrincipalAxes {
            positions,
            masses,
            options,
        } => {
            let positions = coordinates(operation, input.arrays, *positions)?;
            let masses = scalar_values(operation, input.scalars, *masses, positions.len())?;
            pdbiox_geom::principal_axes_with_options(positions, masses, *options)
                .map(GeometryValue::PrincipalAxes)
                .map_err(ExecutionPlanError::Eigen)
        }
        GeometryRequest::Asphericity { positions, options } => {
            let positions = coordinates(operation, input.arrays, *positions)?;
            pdbiox_geom::asphericity_with_options(positions, *options)
                .map(GeometryValue::Asphericity)
                .map_err(ExecutionPlanError::Eigen)
        }
        GeometryRequest::GyrationAxes { positions, options } => {
            let positions = coordinates(operation, input.arrays, *positions)?;
            pdbiox_geom::gyration_axes_with_options(positions, *options)
                .map(GeometryValue::GyrationAxes)
                .map_err(ExecutionPlanError::Eigen)
        }
        GeometryRequest::DistanceMatrix { positions } => {
            pdbiox_geom::distance_matrix(coordinates(operation, input.arrays, *positions)?)
                .map(GeometryValue::DistanceMatrix)
                .map_err(ExecutionPlanError::Matrix)
        }
        GeometryRequest::DistanceMatrixBetween { left, right } => {
            let left = coordinates(operation, input.arrays, *left)?;
            let right = coordinates(operation, input.arrays, *right)?;
            pdbiox_geom::distance_matrix_between(left, right)
                .map(GeometryValue::DistanceMatrix)
                .map_err(ExecutionPlanError::Matrix)
        }
        GeometryRequest::Rmsf { frames } => {
            let frame_input =
                input
                    .frames
                    .get(*frames)
                    .ok_or_else(|| ExecutionPlanError::FrameSlot {
                        operation: operation.into(),
                        slot: *frames,
                    })?;
            let expected = frame_input
                .frame_count
                .checked_mul(frame_input.atom_count)
                .ok_or_else(|| {
                    ExecutionPlanError::InvalidRequest(
                        "RMSF frame dimensions overflow usize".into(),
                    )
                })?;
            if expected != frame_input.positions.len() {
                return Err(ExecutionPlanError::InvalidRequest(
                    "RMSF frame dimensions do not match the borrowed array".into(),
                ));
            }
            let frame_views = if frame_input.atom_count == 0 {
                std::iter::repeat_n(&frame_input.positions[..0], frame_input.frame_count)
                    .collect::<Vec<_>>()
            } else {
                frame_input
                    .positions
                    .chunks_exact(frame_input.atom_count)
                    .collect::<Vec<_>>()
            };
            pdbiox_geom::rmsf(&frame_views)
                .map(GeometryValue::Rmsf)
                .map_err(ExecutionPlanError::Fluctuation)
        }
    }
}

fn coordinates<'a>(
    operation: &str,
    arrays: &'a [CoordinateInput<'a>],
    slot: usize,
) -> Result<&'a [[f32; 3]], ExecutionPlanError> {
    arrays
        .get(slot)
        .map(|input| input.positions)
        .ok_or_else(|| ExecutionPlanError::ArraySlot {
            operation: operation.into(),
            slot,
        })
}

fn scalar_values<'a>(
    operation: &str,
    scalars: &'a [ScalarInput<'a>],
    slot: Option<usize>,
    positions: usize,
) -> Result<&'a [f64], ExecutionPlanError> {
    let Some(slot) = slot else {
        return Ok(&[]);
    };
    let values = scalars
        .get(slot)
        .ok_or_else(|| ExecutionPlanError::ScalarSlot {
            operation: operation.into(),
            slot,
        })?;
    if values.values.len() != positions {
        return Err(ExecutionPlanError::InvalidRequest(
            format!("operation {operation:?} masses must contain exactly one value per position")
                .into(),
        ));
    }
    Ok(values.values)
}
