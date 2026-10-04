//! Running a structure kernel under an explicit policy and wrapping its result.

use crate::analysis_result::PyAnalysis;
use crate::bindings::PyStructure;
use crate::policy::{PyAnalysisPolicy, policy_of};
use crate::table::TableBuilder;
use molframe::analysis::{StructureKernel, analyse_structure};
use pyo3::prelude::*;
use std::fmt::Debug;

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
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<PyAnalysis>
where
    K: StructureKernel,
    K::Error: Debug + std::fmt::Display,
    K::Error: 'static,
    for<'a> molframe::Diagnostic: From<&'a K::Error>,
{
    let structure = structure.inner.clone();
    let analysis = crate::execution::run(py, context, |context| {
        analyse_structure(structure.engine(), policy, kernel, context)
    })?
    .map_err(|error| {
        crate::error::from_diagnostic(&molframe::analysis::governed_diagnostic(&error, |inner| {
            molframe::Diagnostic::from(inner)
        }))
    })?;
    let envelope = PyAnalysis::new(&analysis, None);
    let value = match analysis.into_result() {
        Ok(value) => Some(convert(py, value)?),
        Err(_) => None,
    };
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
    context=None,
))]
pub(crate) fn hydrogen_bonds(
    py: Python<'_>,
    structure: &PyStructure,
    max_distance: f32,
    min_angle: f64,
    backend: &str,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let options = molframe::analysis::HydrogenBondOptions {
        maximum_donor_acceptor_distance: max_distance,
        minimum_angle_degrees: min_angle,
        backend: crate::backend::parse(backend)?,
        periodic: false,
    };
    let kernel = molframe::analysis::hydrogen_bonds_kernel(options);
    run(
        py,
        structure,
        &policy_of(policy),
        &kernel,
        |py, table| {
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
        },
        context,
    )
}

/// Salt bridges between oppositely charged CCD atoms, under a policy.
#[pyfunction]
#[pyo3(signature = (structure, *, max_distance=4.0, backend="auto", policy=None, context=None))]
pub(crate) fn salt_bridges(
    py: Python<'_>,
    structure: &PyStructure,
    max_distance: f32,
    backend: &str,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let kernel =
        molframe::analysis::salt_bridges_kernel(max_distance, crate::backend::parse(backend)?);
    run(
        py,
        structure,
        &policy_of(policy),
        &kernel,
        |py, table| {
            let table = TableBuilder::new(py, table.len())
                .indices("anion", &indices(table.anion()))
                .indices("cation", &indices(table.cation()))
                .single("distance", table.distance())
                .finish();
            Ok(Py::new(py, table)?.into_any())
        },
        context,
    )
}
