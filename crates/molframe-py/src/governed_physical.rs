//! Governed physical analyses of one structure over atom groups: membrane
//! leaflets, radial distributions and coordination.
//!
//! Every binning and every threshold is a required keyword. Atom groups are
//! named by a selection (query text, a compiled query or a selection of the
//! same structure).

use crate::analysis_result::PyAnalysis;
use crate::bindings::{PyQuery, PySelection, PyStructure};
use crate::execution::PyExecutionContext;
use crate::governed::run;
use crate::policy::{PyAnalysisPolicy, policy_of, select_compiled};
use crate::table::TableBuilder;
use molframe::analysis::{CartesianAxis, LeafletOptions, RadialDistributionOptions};
use molframe::engine::core::AtomSelection;
use pyo3::prelude::*;

/// The atoms a selection argument names, in ascending order.
pub(crate) fn atoms_of(
    py: Python<'_>,
    structure: &PyStructure,
    value: &Bound<'_, PyAny>,
    policy: &molframe::AnalysisPolicy,
) -> PyResult<AtomSelection> {
    let count = structure.inner.engine().atom_count();
    let indices = if let Ok(selection) = value.extract::<PyRef<'_, PySelection>>() {
        selection.atom_indices().to_vec()
    } else if let Ok(query) = value.extract::<PyRef<'_, PyQuery>>() {
        select_compiled(py, structure, query.native(), policy)?
            .atom_indices()
            .to_vec()
    } else {
        let source = value.extract::<&str>().map_err(|_| {
            crate::error::type_error("a selection must be query text, a Query or a Selection")
        })?;
        let compiled = molframe::Query::compile(source).map_err(|findings| {
            crate::query_messages::query_error(&molframe::Findings::from(findings), source)
        })?;
        select_compiled(py, structure, &compiled, policy)?
            .atom_indices()
            .to_vec()
    };
    if indices.iter().any(|&atom| atom >= count) {
        return Err(crate::error::value(
            "the selection belongs to a different structure",
        ));
    }
    Ok(AtomSelection::from_sorted(indices))
}

/// Per-atom weights: `"count"` (one each, the default), `"mass"` (atomic weight) or an array.
pub(crate) fn weights_of(
    structure: &PyStructure,
    value: Option<&Bound<'_, PyAny>>,
) -> PyResult<Vec<f64>> {
    let engine = structure.inner.engine();
    let Some(value) = value else {
        return Ok(vec![1.0; engine.atom_count() as usize]);
    };
    if let Ok(name) = value.extract::<&str>() {
        return match name {
            "count" => Ok(vec![1.0; engine.atom_count() as usize]),
            "mass" => {
                let masses = molframe::chemistry::atom_masses(engine);
                let unknown = masses.iter().filter(|mass| mass.is_nan()).count();
                if unknown > 0 {
                    return Err(crate::error::value(format!(
                        "{unknown} atom(s) have an unknown element, so no atomic weight"
                    )));
                }
                Ok(masses)
            }
            other => Err(crate::error::value(format!(
                "weights must be \"count\", \"mass\" or an array, not {other:?}"
            ))),
        };
    }
    value.extract::<Vec<f64>>()
}

pub(crate) fn axis_of(name: &str) -> PyResult<CartesianAxis> {
    name.parse().map_err(crate::error::kernel)
}

/// Connected membrane leaflets of the selected representative sites.
///
/// Two sites are connected when they lie within `connection_distance` ångström;
/// a leaflet is a connected component. Which atom represents a lipid is the
/// caller's choice, so the selection is explicit. Columns: `site`, `leaflet`.
#[pyfunction]
#[pyo3(signature = (
    structure,
    sites,
    *,
    connection_distance,
    backend="auto",
    policy=None,
    context=None,
))]
pub(crate) fn leaflets(
    py: Python<'_>,
    structure: &PyStructure,
    sites: &Bound<'_, PyAny>,
    connection_distance: f32,
    backend: &str,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let policy = policy_of(policy);
    let sites = atoms_of(py, structure, sites, &policy)?;
    let kernel = molframe::analysis::leaflets_kernel(
        &sites,
        LeafletOptions {
            connection_distance,
            backend: crate::backend::parse(backend)?,
        },
    );
    run(
        py,
        structure,
        &policy,
        &kernel,
        |py, found| {
            let mut site = Vec::new();
            let mut leaflet = Vec::new();
            for (number, component) in (0_u32..).zip(&found) {
                for &atom in &component.sites {
                    site.push(atom);
                    leaflet.push(number);
                }
            }
            let table = TableBuilder::new(py, site.len())
                .indices("site", &site)
                .indices("leaflet", &leaflet)
                .finish();
            Ok(Py::new(py, table)?.into_any())
        },
        context,
    )
}

/// The radial distribution function between two selections.
///
/// `volume` (Å³) normalises the shell counts: the ideal population of a shell
/// is `N_pairs * shell_volume / volume`, with `N(N-1)/2` pairs when the two
/// selections are the same atoms. Columns: `lower`, `upper`, `count`,
/// `distribution`.
#[pyfunction]
#[pyo3(signature = (
    structure,
    first,
    second,
    *,
    minimum_distance,
    maximum_distance,
    bins,
    volume,
    backend="auto",
    policy=None,
    context=None,
))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn radial_distribution(
    py: Python<'_>,
    structure: &PyStructure,
    first: &Bound<'_, PyAny>,
    second: &Bound<'_, PyAny>,
    minimum_distance: f32,
    maximum_distance: f32,
    bins: usize,
    volume: f64,
    backend: &str,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let policy = policy_of(policy);
    let left = atoms_of(py, structure, first, &policy)?;
    let right = atoms_of(py, structure, second, &policy)?;
    let kernel = molframe::analysis::radial_distribution_kernel(
        &left,
        &right,
        RadialDistributionOptions {
            minimum_distance,
            maximum_distance,
            bins,
            volume,
            backend: crate::backend::parse(backend)?,
        },
    );
    run(
        py,
        structure,
        &policy,
        &kernel,
        |py, rows| {
            let lower: Vec<f32> = rows.iter().map(|row| row.lower).collect();
            let upper: Vec<f32> = rows.iter().map(|row| row.upper).collect();
            let count: Vec<u64> = rows.iter().map(|row| row.count).collect();
            let distribution: Vec<f64> = rows.iter().map(|row| row.distribution).collect();
            let table = TableBuilder::new(py, rows.len())
                .single("lower", &lower)
                .single("upper", &upper)
                .counts("count", &count)
                .double("distribution", &distribution)
                .finish();
            Ok(Py::new(py, table)?.into_any())
        },
        context,
    )
}

/// Neighbours of every `first` atom among the `second` atoms in a distance
/// shell `[minimum_distance, maximum_distance]`. Columns: `atom`, `count`.
#[pyfunction]
#[pyo3(signature = (
    structure,
    first,
    second,
    *,
    minimum_distance,
    maximum_distance,
    backend="auto",
    policy=None,
    context=None,
))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn coordination_numbers(
    py: Python<'_>,
    structure: &PyStructure,
    first: &Bound<'_, PyAny>,
    second: &Bound<'_, PyAny>,
    minimum_distance: f32,
    maximum_distance: f32,
    backend: &str,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyAnalysis> {
    let policy = policy_of(policy);
    let left = atoms_of(py, structure, first, &policy)?;
    let right = atoms_of(py, structure, second, &policy)?;
    let kernel = molframe::analysis::coordination_numbers_kernel(
        &left,
        &right,
        minimum_distance,
        maximum_distance,
        crate::backend::parse(backend)?,
    );
    let atoms: Vec<u32> = left.into_iter().collect();
    run(
        py,
        structure,
        &policy,
        &kernel,
        |py, counts| {
            let table = TableBuilder::new(py, counts.len())
                .indices("atom", &atoms)
                .indices("count", &counts)
                .finish();
            Ok(Py::new(py, table)?.into_any())
        },
        context,
    )
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(leaflets, module)?)?;
    module.add_function(wrap_pyfunction!(radial_distribution, module)?)?;
    module.add_function(wrap_pyfunction!(coordination_numbers, module)?)
}
