//! Dense benchmark measurements delegated to Rust geometry kernels.

use super::errors::ame_measurement_error;
use crate::geometry::borrowed_coordinates;
use numpy::PyReadonlyArray2;
use pyo3::prelude::*;
use std::collections::BTreeMap;

#[pyclass(name = "MotifBenchMeasurements", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyMotifBenchMeasurements(pub(crate) molframe::fx::MotifBenchMeasurements);

#[pymethods]
impl PyMotifBenchMeasurements {
    #[getter]
    const fn rmsd(&self) -> f64 {
        self.0.rmsd
    }

    #[getter]
    const fn motif_rmsd(&self) -> f64 {
        self.0.motif_rmsd
    }

    fn metrics(&self) -> BTreeMap<String, f64> {
        self.0
            .metrics()
            .into_iter()
            .map(|(name, value)| (name.to_string(), value))
            .collect()
    }
}

#[pyclass(name = "AmeMeasurements", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyAmeMeasurements(pub(crate) molframe::fx::AmeMeasurements);

#[pymethods]
impl PyAmeMeasurements {
    #[getter]
    const fn catalytic_heavy_atom_rmsd(&self) -> f64 {
        self.0.catalytic_heavy_atom_rmsd
    }

    #[getter]
    const fn ligand_backbone_min_distance(&self) -> f64 {
        self.0.ligand_backbone_min_distance
    }

    fn metrics(&self) -> BTreeMap<String, f64> {
        self.0
            .metrics()
            .into_iter()
            .map(|(name, value)| (name.to_string(), value))
            .collect()
    }
}

#[pyfunction]
pub(crate) fn measure_motifbench(
    py: Python<'_>,
    generated_c_alpha: PyReadonlyArray2<'_, f32>,
    predicted_c_alpha: PyReadonlyArray2<'_, f32>,
    reference_motif_backbone: PyReadonlyArray2<'_, f32>,
    predicted_motif_backbone: PyReadonlyArray2<'_, f32>,
) -> PyResult<PyMotifBenchMeasurements> {
    let generated_c_alpha = borrowed_coordinates(&generated_c_alpha)?;
    let predicted_c_alpha = borrowed_coordinates(&predicted_c_alpha)?;
    let reference_motif_backbone = borrowed_coordinates(&reference_motif_backbone)?;
    let predicted_motif_backbone = borrowed_coordinates(&predicted_motif_backbone)?;
    py.detach(|| {
        molframe::fx::measure_motifbench(
            generated_c_alpha,
            predicted_c_alpha,
            reference_motif_backbone,
            predicted_motif_backbone,
        )
        .map(PyMotifBenchMeasurements)
        .map_err(|error| pyo3::exceptions::PyValueError::new_err(format!("{error:?}")))
    })
}

#[pyfunction]
pub(crate) fn measure_ame(
    py: Python<'_>,
    generated_catalytic_backbone: PyReadonlyArray2<'_, f32>,
    predicted_catalytic_backbone: PyReadonlyArray2<'_, f32>,
    generated_catalytic_heavy: PyReadonlyArray2<'_, f32>,
    predicted_catalytic_heavy: PyReadonlyArray2<'_, f32>,
    predicted_ligand: PyReadonlyArray2<'_, f32>,
    predicted_backbone: PyReadonlyArray2<'_, f32>,
) -> PyResult<PyAmeMeasurements> {
    let generated_catalytic_backbone = borrowed_coordinates(&generated_catalytic_backbone)?;
    let predicted_catalytic_backbone = borrowed_coordinates(&predicted_catalytic_backbone)?;
    let generated_catalytic_heavy = borrowed_coordinates(&generated_catalytic_heavy)?;
    let predicted_catalytic_heavy = borrowed_coordinates(&predicted_catalytic_heavy)?;
    let predicted_ligand = borrowed_coordinates(&predicted_ligand)?;
    let predicted_backbone = borrowed_coordinates(&predicted_backbone)?;
    py.detach(|| {
        molframe::fx::measure_ame(
            generated_catalytic_backbone,
            predicted_catalytic_backbone,
            generated_catalytic_heavy,
            predicted_catalytic_heavy,
            predicted_ligand,
            predicted_backbone,
        )
        .map(PyAmeMeasurements)
        .map_err(ame_measurement_error)
    })
}
