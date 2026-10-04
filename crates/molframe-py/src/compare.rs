//! Mechanical adapters for superposition-based structure comparison scores.

use numpy::{PyArray1, PyArray2, PyArrayMethods};
use pyo3::prelude::*;

type Scorer = fn(&[[f32; 3]], &[[f32; 3]]) -> Result<f64, molframe::compare::CompareError>;

fn score(
    py: Python<'_>,
    model: &Bound<'_, PyArray2<f32>>,
    reference: &Bound<'_, PyArray2<f32>>,
    scorer: Scorer,
) -> PyResult<f64> {
    let model = model.readonly();
    let reference = reference.readonly();
    let (model, reference) = (
        crate::bindings::coordinates(&model)?,
        crate::bindings::coordinates(&reference)?,
    );
    py.detach(|| scorer(model, reference))
        .map_err(crate::error::kernel)
}

/// TM-score of `model` fitted onto `reference` (residue-matched coordinates).
#[pyfunction]
fn tm_score(
    py: Python<'_>,
    model: &Bound<'_, PyArray2<f32>>,
    reference: &Bound<'_, PyArray2<f32>>,
) -> PyResult<f64> {
    score(py, model, reference, molframe::compare::tm_score)
}

/// GDT-TS: mean fraction of atoms within 1, 2, 4 and 8 Å after fitting.
#[pyfunction]
fn gdt_ts(
    py: Python<'_>,
    model: &Bound<'_, PyArray2<f32>>,
    reference: &Bound<'_, PyArray2<f32>>,
) -> PyResult<f64> {
    score(py, model, reference, molframe::compare::gdt_ts)
}

/// GDT-HA: mean fraction of atoms within 0.5, 1, 2 and 4 Å after fitting.
#[pyfunction]
fn gdt_ha(
    py: Python<'_>,
    model: &Bound<'_, PyArray2<f32>>,
    reference: &Bound<'_, PyArray2<f32>>,
) -> PyResult<f64> {
    score(py, model, reference, molframe::compare::gdt_ha)
}

/// RMSD of already-aligned coordinates with a non-negative weight per atom.
#[pyfunction]
fn weighted_rmsd(
    py: Python<'_>,
    model: &Bound<'_, PyArray2<f32>>,
    reference: &Bound<'_, PyArray2<f32>>,
    weights: &Bound<'_, PyArray1<f64>>,
) -> PyResult<f64> {
    let model = model.readonly();
    let reference = reference.readonly();
    let weights = weights.readonly();
    let weights = weights.as_slice().map_err(|_| {
        crate::error::value("weights must be C-contiguous; call numpy.ascontiguousarray")
    })?;
    let (model, reference) = (
        crate::bindings::coordinates(&model)?,
        crate::bindings::coordinates(&reference)?,
    );
    py.detach(|| molframe::compare::weighted_rmsd(model, reference, weights))
        .map_err(crate::error::kernel)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    crate::compare_structures::register(module)?;
    module.add_function(wrap_pyfunction!(tm_score, module)?)?;
    module.add_function(wrap_pyfunction!(gdt_ts, module)?)?;
    module.add_function(wrap_pyfunction!(gdt_ha, module)?)?;
    module.add_function(wrap_pyfunction!(weighted_rmsd, module)?)
}
