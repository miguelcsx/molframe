//! Lazy, bounded-memory motif matching: count, first-n, chunked iteration.
//!
//! Matches are residue-index tuples in a compact `(n, k)` `uint32` array whose
//! columns follow the motif's components in sorted-name order. Nothing here
//! builds a Python object per match, and every call counts before it allocates.

use crate::bindings::{PyStructure, findings_error};
use crate::execution::{PyExecutionContext, run};
use crate::policy::{PyAnalysisPolicy, policy_of};
use molframe::motif::{
    EvaluationSpecification, MappingError, MappingPlan, MappingStream,
    read_evaluation_specification,
};
use molframe::{CancellationToken, ExecutionContext};
use numpy::{PyArray2, PyArrayMethods, PyUntypedArrayMethods};
use pyo3::exceptions::PyStopIteration;
use pyo3::prelude::*;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

/// Bytes one match occupies per component column.
const BYTES_PER_INDEX: usize = std::mem::size_of::<u32>();

/// The exception a mapping failure is raised as.
pub(crate) fn mapping_error(error: &MappingError) -> PyErr {
    match error {
        MappingError::Cancelled => crate::error::named("Cancelled", error.to_string()),
        _ => crate::error::failure(error),
    }
}

/// Reads the specification and prepares the mapping plan for `source`.
pub(crate) fn prepare(
    source: &molframe::Structure,
    specification: &std::path::Path,
    components: Option<PathBuf>,
    version: &str,
    policy: &molframe::AnalysisPolicy,
    constrained: bool,
) -> PyResult<(EvaluationSpecification, Arc<MappingPlan>)> {
    let specification =
        read_evaluation_specification(specification).map_err(crate::error::failure)?;
    let provider = components
        .map(|path| {
            molframe::read_component_dictionary(&path, molframe::DictionaryVersion::new(version))
                .map(|(provider, _)| provider)
                .map_err(|findings| findings_error(&findings))
        })
        .transpose()?;
    let plan = MappingPlan::build(
        source.engine(),
        &specification.motif,
        provider
            .as_ref()
            .map(|provider| provider as &dyn molframe::chemistry::ComponentProvider),
        policy,
        constrained,
    )
    .map_err(|error| mapping_error(&error))?;
    Ok((specification, Arc::new(plan)))
}

fn budget_error(what: &str, bytes: usize, context: &ExecutionContext) -> PyErr {
    crate::error::named(
        "MemoryBudgetError",
        format!(
            "{what} needs about {bytes} bytes, over the execution memory budget of {} bytes; \
             raise it with ExecutionContext(memory_budget=<bytes>) passed as context=, \
             lower the number of matches requested, or stream with motif.iter_matches(...) \
             and motif.count(...), which hold only one chunk at a time",
            context.memory_budget().bytes()
        ),
    )
}

/// Matches of a motif as compact residue-index tuples.
#[derive(Debug)]
#[pyclass(name = "MotifMatches", frozen, module = "molframe.motif")]
struct PyMotifMatches {
    components: Vec<String>,
    residues: Py<PyArray2<u32>>,
}

#[pymethods]
impl PyMotifMatches {
    /// Component names, one per column of `residues`.
    #[getter]
    fn components(&self) -> Vec<String> {
        self.components.clone()
    }

    /// A `(matches, components)` `uint32` array of residue indices.
    #[getter]
    fn residues(&self, py: Python<'_>) -> Py<PyArray2<u32>> {
        self.residues.clone_ref(py)
    }

    fn __len__(&self, py: Python<'_>) -> usize {
        self.residues.bind(py).shape()[0]
    }

    fn __repr__(&self, py: Python<'_>) -> String {
        format!(
            "MotifMatches(matches={}, components={:?})",
            self.__len__(py),
            self.components
        )
    }
}

fn names(plan: &MappingPlan) -> Vec<String> {
    plan.components().iter().map(ToString::to_string).collect()
}

fn to_matches(py: Python<'_>, plan: &MappingPlan, flat: Vec<u32>) -> PyResult<PyMotifMatches> {
    let width = plan.width();
    let rows = match flat.len().checked_div(width) {
        Some(rows) => rows,
        None => 0,
    };
    let array = numpy::PyArray1::from_vec(py, flat).reshape([rows, width])?;
    Ok(PyMotifMatches {
        components: names(plan),
        residues: array.unbind(),
    })
}

/// Fills `flat` with up to `limit` tuples from `stream`; `false` when exhausted early.
fn drain(
    plan: &MappingPlan,
    stream: &mut MappingStream,
    limit: usize,
    cancel: &CancellationToken,
) -> Result<Vec<u32>, MappingError> {
    let width = plan.width();
    let mut flat = Vec::new();
    let mut tuple = vec![0u32; width];
    let mut taken = 0usize;
    while taken < limit && stream.next_into(&mut tuple, cancel)? {
        flat.extend_from_slice(&tuple);
        taken += 1;
    }
    Ok(flat)
}

type Outcome<T> = Result<T, PyErr>;

#[allow(clippy::too_many_arguments)]
fn plan_of(
    py: Python<'_>,
    structure: &PyStructure,
    specification: PathBuf,
    components: Option<PathBuf>,
    version: &str,
    constrained: bool,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<Arc<MappingPlan>> {
    let source = structure.inner.clone();
    let policy = policy_of(policy);
    let version = version.to_owned();
    run(py, context, move |_| {
        prepare(
            &source,
            &specification,
            components,
            &version,
            &policy,
            constrained,
        )
        .map(|(_, plan)| plan)
    })?
}

/// Counts the motif's mappings without materialising any.
///
/// With `constrained=True` (the default) only mappings whose distance
/// constraints are satisfied are counted, found with a spatial index instead
/// of a Cartesian product. `constrained=False` counts every structural mapping.
#[pyfunction]
#[pyo3(signature = (structure, specification, *, components=None, version="unversioned", constrained=true, policy=None, context=None))]
#[allow(clippy::too_many_arguments)]
fn count(
    py: Python<'_>,
    structure: &PyStructure,
    specification: PathBuf,
    components: Option<PathBuf>,
    version: &str,
    constrained: bool,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<u64> {
    let plan = plan_of(
        py,
        structure,
        specification,
        components,
        version,
        constrained,
        policy,
        context,
    )?;
    run(py, context, |ctx| plan.count(ctx.cancellation()))?.map_err(|error| mapping_error(&error))
}

/// All matches as one compact array, refused before allocation if over budget.
///
/// The count is taken first; when `n * k * 4` bytes would exceed the context's
/// memory budget, or `limit` mappings, the call raises instead of allocating.
#[pyfunction]
#[pyo3(signature = (structure, specification, *, components=None, version="unversioned", constrained=true, limit=None, policy=None, context=None))]
#[allow(clippy::too_many_arguments)]
fn matches(
    py: Python<'_>,
    structure: &PyStructure,
    specification: PathBuf,
    components: Option<PathBuf>,
    version: &str,
    constrained: bool,
    limit: Option<usize>,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyMotifMatches> {
    let plan = plan_of(
        py,
        structure,
        specification,
        components,
        version,
        constrained,
        policy,
        context,
    )?;
    let flat: Outcome<Vec<u32>> = run(py, context, |ctx| {
        let cancel = ctx_token(ctx);
        let total = plan.count(&cancel).map_err(|error| mapping_error(&error))?;
        let total = match usize::try_from(total) {
            Ok(total) => total,
            Err(_) => usize::MAX,
        };
        if let Some(limit) = limit.filter(|&limit| total > limit) {
            return Err(crate::error::value(format!(
                "the motif has {total} matches, over limit={limit}; raise limit, \
                 or stream with motif.iter_matches(...)"
            )));
        }
        let bytes = total
            .saturating_mul(plan.width())
            .saturating_mul(BYTES_PER_INDEX);
        let _reservation = ctx
            .try_reserve(bytes)
            .map_err(|_| budget_error(&format!("{total} matches"), bytes, ctx))?;
        let mut stream = plan.stream();
        drain(&plan, &mut stream, total, &cancel).map_err(|error| mapping_error(&error))
    })?;
    to_matches(py, &plan, flat?)
}

/// The first `n` matches, in deterministic order, without counting the rest.
#[pyfunction]
#[pyo3(signature = (structure, specification, n, *, components=None, version="unversioned", constrained=true, policy=None, context=None))]
#[allow(clippy::too_many_arguments)]
fn first(
    py: Python<'_>,
    structure: &PyStructure,
    specification: PathBuf,
    n: usize,
    components: Option<PathBuf>,
    version: &str,
    constrained: bool,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyMotifMatches> {
    let plan = plan_of(
        py,
        structure,
        specification,
        components,
        version,
        constrained,
        policy,
        context,
    )?;
    let flat: Outcome<Vec<u32>> = run(py, context, |ctx| {
        let bytes = n
            .saturating_mul(plan.width())
            .saturating_mul(BYTES_PER_INDEX);
        let _reservation = ctx
            .try_reserve(bytes)
            .map_err(|_| budget_error(&format!("first({n})"), bytes, ctx))?;
        let mut stream = plan.stream();
        drain(&plan, &mut stream, n, &ctx_token(ctx)).map_err(|error| mapping_error(&error))
    })?;
    to_matches(py, &plan, flat?)
}

/// The cancellation token an operation's context carries.
fn ctx_token(context: &ExecutionContext) -> CancellationToken {
    context.cancellation().clone()
}

/// Streams matches as `(m, k)` `uint32` arrays of at most `chunk_size` rows.
#[derive(Debug)]
#[pyclass(name = "MatchIterator", frozen, module = "molframe.motif")]
struct PyMatchIterator {
    plan: Arc<MappingPlan>,
    stream: Mutex<MappingStream>,
    chunk_size: usize,
    context: Option<Py<PyExecutionContext>>,
}

#[pymethods]
impl PyMatchIterator {
    /// Component names, one per column of every chunk.
    #[getter]
    fn components(&self) -> Vec<String> {
        names(&self.plan)
    }

    const fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&self, py: Python<'_>) -> PyResult<Py<PyArray2<u32>>> {
        let context = self.context.as_ref().map(Py::get);
        let chunk = self.chunk_size;
        let flat: Outcome<Vec<u32>> = run(py, context, |ctx| {
            let bytes = chunk
                .saturating_mul(self.plan.width())
                .saturating_mul(BYTES_PER_INDEX);
            let _reservation = ctx
                .try_reserve(bytes)
                .map_err(|_| budget_error(&format!("a chunk of {chunk} matches"), bytes, ctx))?;
            let mut stream = self
                .stream
                .lock()
                .map_err(|_| crate::error::internal("match iterator lock poisoned"))?;
            drain(&self.plan, &mut stream, chunk, &ctx_token(ctx))
                .map_err(|error| mapping_error(&error))
        })?;
        let flat = flat?;
        if flat.is_empty() && self.plan.width() > 0 {
            return Err(PyStopIteration::new_err(()));
        }
        if self.plan.width() == 0 {
            // A motif with no components has exactly one empty match.
            return Err(PyStopIteration::new_err(()));
        }
        let rows = flat.len() / self.plan.width();
        Ok(numpy::PyArray1::from_vec(py, flat)
            .reshape([rows, self.plan.width()])?
            .unbind())
    }
}

/// Iterates matches lazily in `(m, k)` `uint32` chunks of at most `chunk_size`.
///
/// Memory stays at one chunk however many matches exist. The chunk is reserved
/// against the context's memory budget on every step.
#[pyfunction]
#[pyo3(signature = (structure, specification, *, chunk_size=65536, components=None, version="unversioned", constrained=true, policy=None, context=None))]
#[allow(clippy::too_many_arguments)]
fn iter_matches(
    py: Python<'_>,
    structure: &PyStructure,
    specification: PathBuf,
    chunk_size: usize,
    components: Option<PathBuf>,
    version: &str,
    constrained: bool,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<Py<PyExecutionContext>>,
) -> PyResult<PyMatchIterator> {
    if chunk_size == 0 {
        return Err(crate::error::value("chunk_size must be positive"));
    }
    let plan = plan_of(
        py,
        structure,
        specification,
        components,
        version,
        constrained,
        policy,
        context.as_ref().map(Py::get),
    )?;
    Ok(PyMatchIterator {
        stream: Mutex::new(plan.stream()),
        plan,
        chunk_size,
        context,
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyMotifMatches>()?;
    module.add_class::<PyMatchIterator>()?;
    module.add_function(wrap_pyfunction!(count, module)?)?;
    module.add_function(wrap_pyfunction!(matches, module)?)?;
    module.add_function(wrap_pyfunction!(first, module)?)?;
    module.add_function(wrap_pyfunction!(iter_matches, module)?)
}
