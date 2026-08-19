//! Lower typed trajectory nodes into facade trajectory requests.

use super::model::PyTrajectoryOperation;
use crate::intrinsic::PyFrameAlignment;
use crate::query::PyAnalysisPolicy;
use numpy::{PyReadonlyArray1, PyReadonlyArray3};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};
use std::collections::HashMap;

pub(crate) fn typed_request(value: &Bound<'_, PyAny>) -> PyResult<Option<PyTrajectoryOperation>> {
    super::classes::typed_classes(value)
}

pub(crate) fn serialized_request(
    config: &Bound<'_, PyDict>,
) -> PyResult<Option<PyTrajectoryOperation>> {
    let Some(tag) = config.get_item("operation")? else {
        return Ok(None);
    };
    let tag: &str = tag.extract()?;
    if !matches!(
        tag,
        "trajectory_rmsd_to_reference"
            | "mean_squared_displacement"
            | "pairwise_fitted_rmsd"
            | "generalized_procrustes_mean"
    ) {
        return Ok(None);
    }
    from_dict(config).map(Some)
}

pub(crate) fn to_native<'py>(
    py: Python<'py>,
    operation: &PyTrajectoryOperation,
    frames: &mut Vec<PyReadonlyArray3<'py, f32>>,
    indices: &mut Vec<PyReadonlyArray1<'py, usize>>,
    frame_slots: &mut HashMap<usize, usize>,
    index_slots: &mut HashMap<usize, usize>,
) -> PyResult<pdbiox::TrajectoryRequest> {
    match operation {
        PyTrajectoryOperation::RmsdToReference {
            frames: value,
            reference,
            alignment,
            policy,
        } => Ok(pdbiox::TrajectoryRequest::RmsdToReference {
            frames: super::super::lowering::retain_frame(py, frames, frame_slots, value)?,
            reference: *reference,
            alignment: (*alignment).into(),
            policy: policy.clone(),
        }),
        PyTrajectoryOperation::MeanSquaredDisplacement {
            frames: value,
            atoms,
            maximum_lag,
            policy,
        } => Ok(pdbiox::TrajectoryRequest::MeanSquaredDisplacement {
            frames: super::super::lowering::retain_frame(py, frames, frame_slots, value)?,
            atoms: atoms
                .as_ref()
                .map(|value| {
                    super::super::lowering::retain_indices(py, indices, index_slots, value)
                })
                .transpose()?,
            maximum_lag: *maximum_lag,
            policy: policy.clone(),
        }),
        PyTrajectoryOperation::PairwiseFittedRmsd {
            frames: value,
            memory_limit,
            policy,
        } => Ok(pdbiox::TrajectoryRequest::PairwiseFittedRmsd {
            frames: super::super::lowering::retain_frame(py, frames, frame_slots, value)?,
            memory_limit: *memory_limit,
            policy: policy.clone(),
        }),
        PyTrajectoryOperation::GeneralizedProcrustesMean {
            frames: value,
            tolerance,
            maximum_iterations,
            policy,
        } => Ok(pdbiox::TrajectoryRequest::GeneralizedProcrustesMean {
            frames: super::super::lowering::retain_frame(py, frames, frame_slots, value)?,
            tolerance: *tolerance,
            maximum_iterations: *maximum_iterations,
            policy: policy.clone(),
        }),
    }
}

pub(crate) fn policy_or_default(policy: Option<PyAnalysisPolicy>) -> pdbiox::AnalysisPolicy {
    policy.map_or_else(pdbiox::AnalysisPolicy::default, |value| value.inner)
}

/// Validates operation controls and borrowed-array contracts at construction.
///
/// The check is constant-time with respect to coordinate data: it validates
/// dtype, shape, contiguity, and scalar controls without walking frames.
pub(crate) fn validate_operation(
    py: Python<'_>,
    operation: &PyTrajectoryOperation,
) -> PyResult<()> {
    let (frame_count, atom_count) = match operation {
        PyTrajectoryOperation::RmsdToReference { frames, .. }
        | PyTrajectoryOperation::MeanSquaredDisplacement { frames, .. }
        | PyTrajectoryOperation::PairwiseFittedRmsd { frames, .. }
        | PyTrajectoryOperation::GeneralizedProcrustesMean { frames, .. } => {
            validate_frames(py, frames)?
        }
    };
    if frame_count == 0 || atom_count == 0 {
        return Err(PyValueError::new_err(
            "trajectory frames must contain at least one frame and one atom",
        ));
    }
    match operation {
        PyTrajectoryOperation::RmsdToReference { reference, .. } => {
            if *reference < frame_count {
                Ok(())
            } else {
                Err(PyValueError::new_err(format!(
                    "reference must index one of {frame_count} trajectory frames"
                )))
            }
        }
        PyTrajectoryOperation::MeanSquaredDisplacement {
            atoms, maximum_lag, ..
        } => {
            if *maximum_lag >= frame_count {
                return Err(PyValueError::new_err(
                    "maximum_lag must be smaller than the frame count",
                ));
            }
            validate_atoms(py, atoms.as_ref())
        }
        PyTrajectoryOperation::PairwiseFittedRmsd { memory_limit, .. } => {
            if *memory_limit == 0 {
                Err(PyValueError::new_err("memory_limit must be positive"))
            } else {
                Ok(())
            }
        }
        PyTrajectoryOperation::GeneralizedProcrustesMean {
            tolerance,
            maximum_iterations,
            ..
        } => {
            if !tolerance.is_finite() || *tolerance <= 0.0 {
                return Err(PyValueError::new_err(
                    "tolerance must be finite and positive",
                ));
            }
            if *maximum_iterations == 0 {
                return Err(PyValueError::new_err("maximum_iterations must be positive"));
            }
            Ok(())
        }
    }
}

fn validate_frames(py: Python<'_>, value: &Py<PyAny>) -> PyResult<(usize, usize)> {
    let array = value.bind(py).extract::<PyReadonlyArray3<'_, f32>>()?;
    let input = crate::intrinsic::borrowed_frame_input(&array)?;
    Ok((input.frame_count, input.atom_count))
}

fn validate_atoms(py: Python<'_>, value: Option<&Py<PyAny>>) -> PyResult<()> {
    let Some(value) = value else {
        return Ok(());
    };
    let indices = value.bind(py).extract::<PyReadonlyArray1<'_, usize>>()?;
    indices.as_slice().map(|_| ()).map_err(|_| {
        PyValueError::new_err(
            "atom indices must be C-contiguous; call numpy.ascontiguousarray explicitly",
        )
    })
}

pub(crate) fn from_dict(config: &Bound<'_, PyDict>) -> PyResult<PyTrajectoryOperation> {
    let operation = super::super::required(config, "operation")?;
    let tag: String = operation.extract()?;
    let frames = super::super::required(config, "frames")?.unbind();
    let policy = match config.get_item("policy")? {
        Some(value) if !value.is_none() => value.extract::<PyAnalysisPolicy>()?.inner,
        _ => pdbiox::AnalysisPolicy::default(),
    };
    let result = match tag.as_str() {
        "trajectory_rmsd_to_reference" => Ok(PyTrajectoryOperation::RmsdToReference {
            frames,
            reference: super::super::required(config, "reference")?.extract()?,
            alignment: match config.get_item("alignment")? {
                Some(value) => value.extract::<PyFrameAlignment>()?,
                None => PyFrameAlignment::Unaligned,
            },
            policy,
        }),
        "mean_squared_displacement" => Ok(PyTrajectoryOperation::MeanSquaredDisplacement {
            frames,
            atoms: config
                .get_item("atoms")?
                .filter(|value| !value.is_none())
                .map(Bound::unbind),
            maximum_lag: super::super::required(config, "maximum_lag")?.extract()?,
            policy,
        }),
        "pairwise_fitted_rmsd" => Ok(PyTrajectoryOperation::PairwiseFittedRmsd {
            frames,
            memory_limit: config
                .get_item("memory_limit")?
                .map_or(Ok(pdbiox::traj::DEFAULT_PAIRWISE_MEMORY_LIMIT), |value| {
                    value.extract()
                })?,
            policy,
        }),
        "generalized_procrustes_mean" => Ok(PyTrajectoryOperation::GeneralizedProcrustesMean {
            frames,
            tolerance: config
                .get_item("tolerance")?
                .map_or(Ok(pdbiox::traj::DEFAULT_PROCRUSTES_TOLERANCE), |value| {
                    value.extract()
                })?,
            maximum_iterations: config.get_item("maximum_iterations")?.map_or(
                Ok(pdbiox::traj::DEFAULT_PROCRUSTES_MAXIMUM_ITERATIONS),
                |value| value.extract(),
            )?,
            policy,
        }),
        _ => Err(pyo3::exceptions::PyValueError::new_err(format!(
            "expected a serialized trajectory operation, found {tag:?}"
        ))),
    }?;
    validate_operation(config.py(), &result)?;
    Ok(result)
}
