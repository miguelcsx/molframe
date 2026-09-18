//! Bindings for the block plan and the reduction policy.

use pyo3::prelude::*;
use pyo3::types::PyModule;

/// How work is divided, independently of how many threads will run it.
#[pyclass(name = "BlockPlan", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyBlockPlan(pub(crate) molframe::core::parallel::BlockPlan);

#[pymethods]
impl PyBlockPlan {
    #[new]
    #[pyo3(signature = (count, block = molframe::core::parallel::DEFAULT_BLOCK_ITEMS))]
    fn new(count: usize, block: usize) -> Self {
        Self(molframe::core::parallel::BlockPlan::new(count, block))
    }

    #[getter]
    fn count(&self) -> usize {
        self.0.count()
    }

    #[getter]
    fn block(&self) -> usize {
        self.0.block()
    }

    #[getter]
    fn blocks(&self) -> usize {
        self.0.blocks()
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// The half-open item range one block covers.
    fn range(&self, index: usize) -> Option<(usize, usize)> {
        self.0.range(index).map(|range| (range.start, range.end))
    }

    /// Every block range, in ascending index order.
    fn ranges(&self) -> Vec<(usize, usize)> {
        self.0
            .ranges()
            .map(|range| (range.start, range.end))
            .collect()
    }

    /// How many workers can be usefully applied.
    fn useful_workers(&self, requested: usize) -> usize {
        self.0.useful_workers(requested)
    }

    fn __repr__(&self) -> String {
        format!(
            "BlockPlan(count={}, block={}, blocks={})",
            self.0.count(),
            self.0.block(),
            self.0.blocks()
        )
    }
}

/// Whether a parallel reduction may reassociate.
#[pyclass(name = "ReductionPolicy", eq, frozen, from_py_object)]
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum PyReductionPolicy {
    /// Partials combine in ascending block index, whatever order they finish in.
    Deterministic,
    /// Partials combine as they complete, permitting work stealing.
    Fast,
}

#[pymethods]
impl PyReductionPolicy {
    /// Returns true when combination order is fixed by block index.
    fn is_deterministic(&self) -> bool {
        molframe::core::parallel::ReductionPolicy::from(*self).is_deterministic()
    }

    /// The name recorded in provenance.
    #[getter]
    fn name(&self) -> &'static str {
        molframe::core::parallel::ReductionPolicy::from(*self).name()
    }

    fn __repr__(&self) -> String {
        format!("ReductionPolicy.{}", self.name())
    }
}

impl From<PyReductionPolicy> for molframe::core::parallel::ReductionPolicy {
    fn from(value: PyReductionPolicy) -> Self {
        match value {
            PyReductionPolicy::Deterministic => Self::Deterministic,
            PyReductionPolicy::Fast => Self::Fast,
        }
    }
}

impl From<molframe::core::parallel::ReductionPolicy> for PyReductionPolicy {
    fn from(value: molframe::core::parallel::ReductionPolicy) -> Self {
        match value {
            molframe::core::parallel::ReductionPolicy::Deterministic => Self::Deterministic,
            molframe::core::parallel::ReductionPolicy::Fast => Self::Fast,
        }
    }
}

/// Bounded native scheduling; Python callback heap use is measured separately.
#[pyfunction]
fn try_for_each_block_in(
    py: Python<'_>,
    plan: PyBlockPlan,
    context: &super::execution::PyExecutionContext,
    bytes_per_block: usize,
    block: Py<PyAny>,
    consume: Py<PyAny>,
) -> PyResult<()> {
    let context = context.native();
    py.detach(move || {
        molframe::core::parallel::try_for_each_block_in(
            plan.0,
            &context,
            bytes_per_block,
            |index, range| Python::attach(|py| block.call1(py, (index, (range.start, range.end)))),
            |value| Python::attach(|py| consume.call1(py, (value,))).map(|_| ()),
        )
    })
    .map_err(|error| match error {
        molframe::core::parallel::BlockExecutionError::Operation(error) => error,
        other => pyo3::exceptions::PyRuntimeError::new_err(other.to_string()),
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add(
        "DEFAULT_BLOCK_ITEMS",
        molframe::core::parallel::DEFAULT_BLOCK_ITEMS,
    )?;
    module.add_class::<PyBlockPlan>()?;
    module.add_function(wrap_pyfunction!(try_for_each_block_in, module)?)?;
    module.add_class::<PyReductionPolicy>()?;
    Ok(())
}
