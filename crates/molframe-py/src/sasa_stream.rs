//! Streamed solvent-accessible surface area of a trajectory.
//!
//! Rust chooses the atoms and radii from the topology under the policy's radii
//! rule; this layer converts the arguments and hands the arrays back.

use crate::analysis_result::PyAnalysis;
use crate::bindings::PyStructure;
use crate::governed_field::sasa_options;
use crate::policy::{PyAnalysisPolicy, policy_of};
use crate::trajectory_analyses::PyFrameSeries;
use crate::trajectory_stream::Source;
use molframe::analysis::SasaStreamPlan;
use pyo3::prelude::*;
use pyo3::types::PySlice;
use std::path::PathBuf;

/// Total solvent-accessible area of every streamed frame, in square angstroms.
///
/// `topology` names the atoms and their radii; Rust applies the same rules as
/// `analysis.sasa` (the policy's `vdw_radii`, `E7102`/`E7103`, hydrogens and
/// other conformers left out) and reads only the measured atoms from the
/// trajectory. Only Shrake-Rupley is streamed.
#[pyfunction]
#[pyo3(signature = (
    path,
    topology,
    *,
    probe=None,
    method=None,
    n_points=None,
    hetero=None,
    metals=None,
    radii=None,
    policy=None,
    format=None,
    frames=None,
    memory_limit=None,
    context=None,
))]
#[pyo3(name = "sasa_stream")]
pub(crate) fn sasa_stream_py(
    py: Python<'_>,
    path: PathBuf,
    topology: &PyStructure,
    probe: Option<f32>,
    method: Option<&str>,
    n_points: Option<u16>,
    hetero: Option<&str>,
    metals: Option<&str>,
    radii: Option<&str>,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    format: Option<&str>,
    frames: Option<&Bound<'_, PySlice>>,
    memory_limit: Option<usize>,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let (selection, source) = crate::tables::radii_under(radii, policy.as_deref())?;
    let options = sasa_options(method, n_points, probe, hetero, metals, &source)?;
    let policy = policy_of(policy);
    let plan = SasaStreamPlan::new(topology.inner.engine(), &options, &selection, &policy)
        .map_err(crate::error::kernel)?;
    let trajectory = Source::new(path, format, memory_limit)?;
    let atoms = trajectory.atom_count()?;
    let (reader, bytes) = trajectory.open(frames, Some(plan.atoms().to_vec()))?;
    let analysis = crate::execution::run(py, context, move |context| {
        let mut reader = reader;
        plan.run(&mut reader, atoms, context, bytes)
    })?
    .map_err(crate::error::kernel)?;
    let envelope = PyAnalysis::new(&analysis, None);
    let value = match analysis.into_result() {
        Ok(series) => {
            let series = PyFrameSeries::from_areas(py, series.frames, series.times, series.areas);
            Some(Py::new(py, series)?.into_any())
        }
        Err(_) => None,
    };
    Ok(envelope.with_value(value))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(sasa_stream_py, module)?)
}
