//! Running a structure kernel under an explicit policy and wrapping its result.

use crate::analysis_result::PyAnalysis;
use crate::bindings::PyStructure;
use crate::policy::PyAnalysisPolicy;
use crate::table::TableBuilder;
use molframe::analysis::{StructureKernel, analyse_structure};
use pyo3::{exceptions::PyValueError, prelude::*};
use std::fmt::Debug;

/// The policy a call names, or the default profile.
pub(crate) fn policy_of(policy: Option<PyRef<'_, PyAnalysisPolicy>>) -> molframe::AnalysisPolicy {
    policy.map_or_else(molframe::AnalysisPolicy::default, |policy| policy.0.clone())
}

/// Runs `kernel` on the first model under `policy` and converts the value.
///
/// Alternate conformations are resolved once under the policy, so the result's
/// coverage, status and assumptions describe what was actually analysed.
pub(crate) fn run<K>(
    py: Python<'_>,
    structure: &PyStructure,
    policy: &molframe::AnalysisPolicy,
    kernel: &K,
    convert: impl FnOnce(Python<'_>, K::Output) -> PyResult<Py<PyAny>>,
) -> PyResult<PyAnalysis>
where
    K: StructureKernel,
    K::Error: Debug,
{
    let structure = structure.inner.clone();
    let analysis = py
        .detach(|| {
            analyse_structure(
                structure.engine(),
                policy,
                kernel,
                &molframe::ExecutionContext::default(),
            )
        })
        .map_err(|error| PyValueError::new_err(format!("{error:?}")))?;
    let envelope = PyAnalysis::new(&analysis, py.None());
    let value = convert(py, analysis.value)?;
    Ok(envelope.with_value(value))
}

fn indices(column: &[molframe::AtomIndex]) -> Vec<u32> {
    column.iter().map(|index| index.get()).collect()
}

/// Hydrogen bonds from CCD donor and acceptor annotations, under a policy.
#[pyfunction]
#[pyo3(signature = (
    structure,
    *,
    max_distance=3.5,
    min_angle=120.0,
    backend="auto",
    policy=None,
))]
pub(crate) fn hydrogen_bonds(
    py: Python<'_>,
    structure: &PyStructure,
    max_distance: f32,
    min_angle: f64,
    backend: &str,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
) -> PyResult<PyAnalysis> {
    let options = molframe::analysis::HydrogenBondOptions {
        maximum_donor_acceptor_distance: max_distance,
        minimum_angle_degrees: min_angle,
        backend: crate::backend::parse(backend)?,
        periodic: false,
    };
    let kernel = molframe::analysis::hydrogen_bonds_kernel(options);
    run(py, structure, &policy_of(policy), &kernel, |py, table| {
        let rows = table.len();
        let table = TableBuilder::new(py, rows)
            .indices("donor", &indices(table.donor()))
            .indices("hydrogen", &indices(table.hydrogen()))
            .indices("acceptor", &indices(table.acceptor()))
            .single("donor_acceptor_distance", table.donor_acceptor_distance())
            .single(
                "hydrogen_acceptor_distance",
                table.hydrogen_acceptor_distance(),
            )
            .double("angle_degrees", table.angle_degrees())
            .finish();
        Ok(Py::new(py, table)?.into_any())
    })
}

/// Salt bridges between oppositely charged CCD atoms, under a policy.
#[pyfunction]
#[pyo3(signature = (structure, *, max_distance=4.0, backend="auto", policy=None))]
pub(crate) fn salt_bridges(
    py: Python<'_>,
    structure: &PyStructure,
    max_distance: f32,
    backend: &str,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
) -> PyResult<PyAnalysis> {
    let kernel =
        molframe::analysis::salt_bridges_kernel(max_distance, crate::backend::parse(backend)?);
    run(py, structure, &policy_of(policy), &kernel, |py, table| {
        let table = TableBuilder::new(py, table.len())
            .indices("anion", &indices(table.anion()))
            .indices("cation", &indices(table.cation()))
            .single("distance", table.distance())
            .finish();
        Ok(Py::new(py, table)?.into_any())
    })
}
