//! Mechanical adapters for solvent-accessible area and cavity detection.

use numpy::{PyArray1, PyArray2, PyArrayMethods, PyReadonlyArray1, ToPyArray};
use pyo3::prelude::*;

fn radii<'a>(array: &'a PyReadonlyArray1<'_, f32>) -> PyResult<&'a [f32]> {
    array.as_slice().map_err(|_| {
        crate::error::value("radii must be C-contiguous; call numpy.ascontiguousarray")
    })
}

/// Per-atom solvent-accessible area in square ångström (Shrake–Rupley).
#[pyfunction]
#[pyo3(signature = (coordinates, radii, *, probe=1.4, points=960, context=None))]
fn sasa<'py>(
    py: Python<'py>,
    coordinates: &Bound<'py, PyArray2<f32>>,
    radii: &Bound<'py, PyArray1<f32>>,
    probe: f32,
    points: u16,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let positions = coordinates.readonly();
    let radii = radii.readonly();
    let (positions, radii) = (
        crate::bindings::coordinates(&positions)?,
        self::radii(&radii)?,
    );
    let areas = crate::execution::run(py, context, |context| {
        molframe::surface::shrake_rupley(positions, radii, probe, points, context)
    })?
    .map_err(crate::error::kernel)?;
    Ok(areas.to_pyarray(py))
}

/// Per-atom solvent-accessible area in square ångström (Lee–Richards slices).
#[pyfunction]
#[pyo3(signature = (coordinates, radii, *, probe=1.4, slices=20, context=None))]
fn lee_richards<'py>(
    py: Python<'py>,
    coordinates: &Bound<'py, PyArray2<f32>>,
    radii: &Bound<'py, PyArray1<f32>>,
    probe: f32,
    slices: u16,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let positions = coordinates.readonly();
    let radii = radii.readonly();
    let (positions, radii) = (
        crate::bindings::coordinates(&positions)?,
        self::radii(&radii)?,
    );
    let areas = crate::execution::run(py, context, |context| {
        molframe::surface::lee_richards(positions, radii, probe, slices, context)
    })?
    .map_err(crate::error::kernel)?;
    Ok(areas.to_pyarray(py))
}

/// Enclosed cavities, largest first, as `(volume, (x, y, z), cells)` rows.
#[pyfunction]
#[pyo3(signature = (coordinates, radii, *, probe=1.4, resolution=0.5))]
fn cavities(
    py: Python<'_>,
    coordinates: &Bound<'_, PyArray2<f32>>,
    radii: &Bound<'_, PyArray1<f32>>,
    probe: f32,
    resolution: f32,
) -> PyResult<Vec<(f64, [f32; 3], usize)>> {
    let positions = coordinates.readonly();
    let radii = radii.readonly();
    let (positions, radii) = (
        crate::bindings::coordinates(&positions)?,
        self::radii(&radii)?,
    );
    let found = py
        .detach(|| molframe::surface::cavities(positions, radii, probe, resolution))
        .map_err(crate::error::kernel)?;
    Ok(found
        .into_iter()
        .map(|cavity| (cavity.volume, cavity.representative, cavity.cells))
        .collect())
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(sasa, module)?)?;
    module.add_function(wrap_pyfunction!(lee_richards, module)?)?;
    module.add_function(wrap_pyfunction!(cavities, module)?)
}
