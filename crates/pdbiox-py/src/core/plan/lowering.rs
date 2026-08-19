//! Lower Python operation nodes into typed native facade requests.

use super::{Operation, PyAny, PyReadonlyArray2, geometry, spatial, trajectory};
use crate::geometry::borrowed_coordinates;
use pyo3::prelude::*;
use std::collections::HashMap;

pub(super) fn lower_operation<'py>(
    native: &mut pdbiox::Plan,
    name: &str,
    operation: &Operation,
    py: Python<'py>,
    arrays: &mut Vec<PyReadonlyArray2<'py, f32>>,
    scalar_arrays: &mut Vec<numpy::PyReadonlyArray1<'py, f64>>,
    float_arrays: &mut Vec<numpy::PyReadonlyArray1<'py, f32>>,
    mask_arrays: &mut Vec<numpy::PyReadonlyArray1<'py, bool>>,
    frame_arrays: &mut Vec<numpy::PyReadonlyArray3<'py, f32>>,
    index_arrays: &mut Vec<numpy::PyReadonlyArray1<'py, usize>>,
    coordinate_slots: &mut HashMap<usize, usize>,
    scalar_slots: &mut HashMap<usize, usize>,
    float_slots: &mut HashMap<usize, usize>,
    mask_slots: &mut HashMap<usize, usize>,
    frame_slots: &mut HashMap<usize, usize>,
    index_slots: &mut HashMap<usize, usize>,
) -> PyResult<()> {
    match operation {
        Operation::Contacts(operation) => lower_contacts(native, name, operation)?,
        Operation::Rmsd(operation) => {
            lower_rmsd(native, name, operation, py, arrays, coordinate_slots)?;
        }
        Operation::Comparison(operation) => {
            lower_comparison(native, name, operation, py, arrays, coordinate_slots)?;
        }
        Operation::BasePairs(operation) => lower_structure(native, name, &operation.request)?,
        Operation::Structure(operation) => lower_structure(native, name, operation)?,
        Operation::Selection(operation) => native
            .add_selection(name.to_owned(), operation.request.clone())
            .map_err(plan_error)?,
        Operation::BondInference(operation) => native
            .add_bond_inference(name.to_owned(), operation.options.0)
            .map_err(plan_error)?,
        Operation::Geometry(operation) => lower_geometry(
            native,
            name,
            operation,
            py,
            arrays,
            scalar_arrays,
            frame_arrays,
            coordinate_slots,
            frame_slots,
        )?,
        Operation::Spatial(operation) => {
            lower_spatial(native, name, operation, py, arrays, coordinate_slots)?;
        }
        Operation::Physical(operation) => lower_physical(
            native,
            name,
            operation,
            py,
            scalar_arrays,
            float_arrays,
            scalar_slots,
            float_slots,
        )?,
        Operation::Surface(operation) => lower_surface(
            native,
            name,
            operation,
            py,
            arrays,
            float_arrays,
            mask_arrays,
            coordinate_slots,
            float_slots,
            mask_slots,
        )?,
        Operation::Trajectory(operation) => lower_trajectory(
            native,
            name,
            operation,
            py,
            frame_arrays,
            index_arrays,
            frame_slots,
            index_slots,
        )?,
    }
    Ok(())
}

fn lower_contacts(
    native: &mut pdbiox::Plan,
    name: &str,
    operation: &super::PyContacts,
) -> PyResult<()> {
    let request = pdbiox::ContactsRequest::new(
        operation.left.clone(),
        operation.right.clone(),
        operation.cutoff,
        operation.backend.into(),
        operation.policy.inner.clone(),
    )
    .map_err(plan_error)?;
    native
        .add_contacts(name.to_owned(), request)
        .map_err(plan_error)
}

fn lower_rmsd<'py>(
    native: &mut pdbiox::Plan,
    name: &str,
    operation: &super::PyRmsd,
    py: Python<'py>,
    arrays: &mut Vec<PyReadonlyArray2<'py, f32>>,
    coordinate_slots: &mut HashMap<usize, usize>,
) -> PyResult<()> {
    let (mobile, reference) = retain_coordinate_pair(
        py,
        arrays,
        coordinate_slots,
        &operation.mobile,
        &operation.reference,
    )?;
    native
        .add_rmsd(name.to_owned(), pdbiox::RmsdRequest::new(mobile, reference))
        .map_err(plan_error)
}

fn lower_comparison<'py>(
    native: &mut pdbiox::Plan,
    name: &str,
    operation: &super::PyComparison,
    py: Python<'py>,
    arrays: &mut Vec<PyReadonlyArray2<'py, f32>>,
    coordinate_slots: &mut HashMap<usize, usize>,
) -> PyResult<()> {
    let (mobile, reference) = retain_coordinate_pair(
        py,
        arrays,
        coordinate_slots,
        &operation.mobile,
        &operation.reference,
    )?;
    native
        .add(
            name.to_owned(),
            pdbiox::PlanOperation::Comparison(pdbiox::ComparisonRequest::new(
                mobile,
                reference,
                operation.metric.native(),
            )),
        )
        .map_err(plan_error)
}

fn lower_structure(
    native: &mut pdbiox::Plan,
    name: &str,
    request: &pdbiox::StructureRequest,
) -> PyResult<()> {
    native
        .add(
            name.to_owned(),
            pdbiox::PlanOperation::Structure(Box::new(request.clone())),
        )
        .map_err(plan_error)
}

fn lower_physical<'py>(
    native: &mut pdbiox::Plan,
    name: &str,
    operation: &super::physical::PyPhysicalOperation,
    py: Python<'py>,
    scalar_arrays: &mut Vec<numpy::PyReadonlyArray1<'py, f64>>,
    float_arrays: &mut Vec<numpy::PyReadonlyArray1<'py, f32>>,
    scalar_slots: &mut HashMap<usize, usize>,
    float_slots: &mut HashMap<usize, usize>,
) -> PyResult<()> {
    let request = super::physical::to_native(
        py,
        operation,
        scalar_arrays,
        float_arrays,
        scalar_slots,
        float_slots,
    )?;
    native
        .add(
            name.to_owned(),
            pdbiox::PlanOperation::Physical(Box::new(request)),
        )
        .map_err(plan_error)
}

fn lower_surface<'py>(
    native: &mut pdbiox::Plan,
    name: &str,
    operation: &super::surface::PySurfaceOperation,
    py: Python<'py>,
    arrays: &mut Vec<PyReadonlyArray2<'py, f32>>,
    float_arrays: &mut Vec<numpy::PyReadonlyArray1<'py, f32>>,
    mask_arrays: &mut Vec<numpy::PyReadonlyArray1<'py, bool>>,
    coordinate_slots: &mut HashMap<usize, usize>,
    float_slots: &mut HashMap<usize, usize>,
    mask_slots: &mut HashMap<usize, usize>,
) -> PyResult<()> {
    let request = super::surface::to_native(
        py,
        operation,
        arrays,
        float_arrays,
        mask_arrays,
        coordinate_slots,
        float_slots,
        mask_slots,
    )?;
    native
        .add(
            name.to_owned(),
            pdbiox::PlanOperation::Surface(Box::new(request)),
        )
        .map_err(plan_error)
}

fn lower_geometry<'py>(
    native: &mut pdbiox::Plan,
    name: &str,
    operation: &geometry::PyGeometryOperation,
    py: Python<'py>,
    arrays: &mut Vec<PyReadonlyArray2<'py, f32>>,
    scalar_arrays: &mut Vec<numpy::PyReadonlyArray1<'py, f64>>,
    frame_arrays: &mut Vec<numpy::PyReadonlyArray3<'py, f32>>,
    coordinate_slots: &mut HashMap<usize, usize>,
    frame_slots: &mut HashMap<usize, usize>,
) -> PyResult<()> {
    let request = geometry::to_native(
        py,
        operation,
        arrays,
        scalar_arrays,
        frame_arrays,
        coordinate_slots,
        frame_slots,
    )?;
    native
        .add(
            name.to_owned(),
            pdbiox::PlanOperation::Geometry(Box::new(request)),
        )
        .map_err(plan_error)
}

fn lower_trajectory<'py>(
    native: &mut pdbiox::Plan,
    name: &str,
    operation: &trajectory::PyTrajectoryOperation,
    py: Python<'py>,
    frame_arrays: &mut Vec<numpy::PyReadonlyArray3<'py, f32>>,
    index_arrays: &mut Vec<numpy::PyReadonlyArray1<'py, usize>>,
    frame_slots: &mut HashMap<usize, usize>,
    index_slots: &mut HashMap<usize, usize>,
) -> PyResult<()> {
    let request = trajectory::to_native(
        py,
        operation,
        frame_arrays,
        index_arrays,
        frame_slots,
        index_slots,
    )?;
    native
        .add_trajectory(name.to_owned(), request)
        .map_err(plan_error)
}

fn lower_spatial<'py>(
    native: &mut pdbiox::Plan,
    name: &str,
    operation: &spatial::PySpatialOperation,
    py: Python<'py>,
    arrays: &mut Vec<PyReadonlyArray2<'py, f32>>,
    coordinate_slots: &mut HashMap<usize, usize>,
) -> PyResult<()> {
    let request = spatial::to_native(py, operation, arrays, coordinate_slots)?;
    native
        .add(
            name.to_owned(),
            pdbiox::PlanOperation::Spatial(Box::new(request)),
        )
        .map_err(plan_error)
}

fn retain_coordinate_pair<'py>(
    py: Python<'py>,
    arrays: &mut Vec<PyReadonlyArray2<'py, f32>>,
    coordinate_slots: &mut HashMap<usize, usize>,
    mobile: &Py<PyAny>,
    reference: &Py<PyAny>,
) -> PyResult<(usize, usize)> {
    let mobile_slot = retain_coordinate(py, arrays, coordinate_slots, mobile)?;
    let reference_slot = retain_coordinate(py, arrays, coordinate_slots, reference)?;
    Ok((mobile_slot, reference_slot))
}

pub(super) fn retain_coordinate<'py>(
    py: Python<'py>,
    arrays: &mut Vec<PyReadonlyArray2<'py, f32>>,
    coordinate_slots: &mut HashMap<usize, usize>,
    value: &Py<PyAny>,
) -> PyResult<usize> {
    let key = value.bind(py).as_ptr() as usize;
    if let Some(slot) = coordinate_slots.get(&key).copied() {
        return Ok(slot);
    }
    let array = value.bind(py).extract::<PyReadonlyArray2<'py, f32>>()?;
    borrowed_coordinates(&array)?;
    let slot = arrays.len();
    arrays.push(array);
    coordinate_slots.insert(key, slot);
    Ok(slot)
}

pub(super) fn retain_frame<'py>(
    py: Python<'py>,
    frames: &mut Vec<numpy::PyReadonlyArray3<'py, f32>>,
    frame_slots: &mut HashMap<usize, usize>,
    value: &Py<PyAny>,
) -> PyResult<usize> {
    let key = value.bind(py).as_ptr() as usize;
    if let Some(slot) = frame_slots.get(&key).copied() {
        return Ok(slot);
    }
    let array = value
        .bind(py)
        .extract::<numpy::PyReadonlyArray3<'py, f32>>()?;
    crate::intrinsic::borrowed_frame_input(&array)?;
    let slot = frames.len();
    frames.push(array);
    frame_slots.insert(key, slot);
    Ok(slot)
}

pub(super) fn retain_indices<'py>(
    py: Python<'py>,
    indices: &mut Vec<numpy::PyReadonlyArray1<'py, usize>>,
    index_slots: &mut HashMap<usize, usize>,
    value: &Py<PyAny>,
) -> PyResult<usize> {
    let key = value.bind(py).as_ptr() as usize;
    if let Some(slot) = index_slots.get(&key).copied() {
        return Ok(slot);
    }
    let array = value
        .bind(py)
        .extract::<numpy::PyReadonlyArray1<'py, usize>>()?;
    array.as_slice().map_err(|_| {
        pyo3::exceptions::PyValueError::new_err(
            "atom indices must be C-contiguous; call numpy.ascontiguousarray explicitly",
        )
    })?;
    let slot = indices.len();
    indices.push(array);
    index_slots.insert(key, slot);
    Ok(slot)
}

pub(super) fn plan_error(error: pdbiox::ExecutionPlanError) -> PyErr {
    pyo3::exceptions::PyValueError::new_err(error.to_string())
}
