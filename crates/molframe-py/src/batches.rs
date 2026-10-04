//! Bounded batches of structure rows, for files too large to hold at once.

use crate::execution::PyExecutionContext;
use molframe::{
    Backpressure, BatchDemand, BatchSource, ExecutionContext, Presence, StructureBatchReader,
};
use numpy::{PyArray1, PyArray2, PyArrayMethods, ToPyArray};
use pyo3::prelude::*;
use pyo3::types::PyList;
use std::path::PathBuf;
use std::sync::Mutex;

/// A recorded value, or `NaN` where the file records none.
fn recorded(values: &[f32], presence: &[Presence]) -> Vec<f32> {
    values
        .iter()
        .zip(presence)
        .map(|(value, presence)| {
            if presence.is_present() {
                *value
            } else {
                f32::NAN
            }
        })
        .collect()
}

/// One bounded run of atom rows, in file order.
///
/// Columns are arrays of one entry per row. Names (`chains`, `components`,
/// `atom_names`) are lists of text, the costly columns; `occupancies` and
/// `b_factors` are `NaN` where the file records no value, never a made-up one.
#[derive(Debug)]
#[pyclass(
    name = "StructureBatch",
    frozen,
    skip_from_py_object,
    module = "molframe"
)]
pub(crate) struct PyStructureBatch {
    models: Py<PyArray1<i32>>,
    sequences: Py<PyArray1<i32>>,
    elements: Py<PyArray1<u8>>,
    positions: Py<PyArray2<f32>>,
    occupancies: Py<PyArray1<f32>>,
    b_factors: Py<PyArray1<f32>>,
    atom_site_ids: Py<PyArray1<u32>>,
    heterogens: Py<PyArray1<u8>>,
    chains: Vec<String>,
    components: Vec<String>,
    atom_names: Vec<String>,
    diagnostics: Py<PyList>,
}

macro_rules! array_getters {
    ($($name:ident: $type:ty),+ $(,)?) => {
        #[pymethods]
        impl PyStructureBatch {
            $(
                #[getter]
                fn $name(&self, py: Python<'_>) -> Py<$type> {
                    self.$name.clone_ref(py)
                }
            )+

            /// Chain identifier of each row.
            #[getter]
            fn chains(&self) -> Vec<String> {
                self.chains.clone()
            }

            /// Residue component name of each row.
            #[getter]
            fn components(&self) -> Vec<String> {
                self.components.clone()
            }

            /// Atom name of each row.
            #[getter]
            fn atom_names(&self) -> Vec<String> {
                self.atom_names.clone()
            }

            /// What was wrong with the rows of this batch.
            #[getter]
            fn diagnostics(&self, py: Python<'_>) -> Py<PyList> {
                self.diagnostics.clone_ref(py)
            }

            fn __len__(&self, py: Python<'_>) -> usize {
                numpy::PyUntypedArrayMethods::len(self.elements.bind(py))
            }
        }
    };
}

array_getters!(
    models: PyArray1<i32>,
    sequences: PyArray1<i32>,
    elements: PyArray1<u8>,
    positions: PyArray2<f32>,
    occupancies: PyArray1<f32>,
    b_factors: PyArray1<f32>,
    atom_site_ids: PyArray1<u32>,
    heterogens: PyArray1<u8>,
);

/// The columns of one batch, copied out so the lease can be released.
struct Columns {
    models: Vec<i32>,
    sequences: Vec<i32>,
    elements: Vec<u8>,
    positions: Vec<[f32; 3]>,
    occupancies: Vec<f32>,
    b_factors: Vec<f32>,
    atom_site_ids: Vec<u32>,
    heterogens: Vec<u8>,
    chains: Vec<String>,
    components: Vec<String>,
    atom_names: Vec<String>,
    diagnostics: Vec<molframe::Diagnostic>,
}

/// Iterates the batches of a structural file within the bounds it was opened with.
#[derive(Debug)]
#[pyclass(
    name = "StructureBatches",
    frozen,
    skip_from_py_object,
    module = "molframe"
)]
pub(crate) struct PyStructureBatches {
    source: Mutex<StructureBatchReader>,
    context: ExecutionContext,
    demand: BatchDemand,
}

#[pymethods]
impl PyStructureBatches {
    const fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&self, py: Python<'_>) -> PyResult<Option<PyStructureBatch>> {
        let pulled = py.detach(|| {
            let Ok(mut source) = self.source.lock() else {
                return Err(None);
            };
            match source.next_batch(self.demand, &self.context) {
                Ok(Backpressure::Finished) => Ok(None),
                Ok(Backpressure::Pending) => Err(Some(
                    molframe::Diagnostic::new(molframe::Code::E7001)
                        .with_message("the next batch cannot fit the remaining execution budget"),
                )),
                Ok(Backpressure::Ready(lease)) => {
                    let batch = lease.batch();
                    let text = |ids: &[molframe::SymbolId]| {
                        ids.iter()
                            .map(|id| batch.dictionary().resolve(*id).map(str::to_owned))
                            .map(|name| match name {
                                Some(name) => name,
                                None => String::new(),
                            })
                            .collect::<Vec<_>>()
                    };
                    let (occupancy, occupancy_presence) = batch.occupancies();
                    let (b_factor, b_factor_presence) = batch.b_factors();
                    Ok(Some(Columns {
                        models: batch.models().to_vec(),
                        sequences: batch.sequences().to_vec(),
                        elements: batch.elements().to_vec(),
                        positions: batch.positions().to_vec(),
                        occupancies: recorded(occupancy, occupancy_presence),
                        b_factors: recorded(b_factor, b_factor_presence),
                        atom_site_ids: batch.atom_site_ids().to_vec(),
                        heterogens: batch.heterogens().to_vec(),
                        chains: text(batch.chains()),
                        components: text(batch.components()),
                        atom_names: text(batch.atoms()),
                        diagnostics: batch.diagnostics().to_vec(),
                    }))
                }
                Err(error) => Err(Some(molframe::Diagnostic::from(&error))),
            }
        });
        let columns = match pulled {
            Ok(Some(columns)) => columns,
            Ok(None) => return Ok(None),
            Err(Some(diagnostic)) => return Err(crate::error::from_diagnostic(&diagnostic)),
            Err(None) => {
                return Err(crate::error::internal(
                    "a batch reader was poisoned by a failed read",
                ));
            }
        };
        let rows = columns.positions.len();
        let flat: Vec<f32> = columns.positions.iter().flatten().copied().collect();
        Ok(Some(PyStructureBatch {
            models: columns.models.to_pyarray(py).unbind(),
            sequences: columns.sequences.to_pyarray(py).unbind(),
            elements: columns.elements.to_pyarray(py).unbind(),
            positions: PyArray1::from_vec(py, flat).reshape([rows, 3])?.unbind(),
            occupancies: columns.occupancies.to_pyarray(py).unbind(),
            b_factors: columns.b_factors.to_pyarray(py).unbind(),
            atom_site_ids: columns.atom_site_ids.to_pyarray(py).unbind(),
            heterogens: columns.heterogens.to_pyarray(py).unbind(),
            chains: columns.chains,
            components: columns.components,
            atom_names: columns.atom_names,
            diagnostics: crate::error::diagnostic_list(py, &columns.diagnostics)?.unbind(),
        }))
    }
}

/// Opens a structural file for reading a bounded run of atom rows at a time.
///
/// Memory stays bounded by `rows` and `bytes` per batch and by the `context`
/// budget, whatever the file's size. Formats without a bounded reader are
/// refused rather than read whole.
#[pyfunction]
#[pyo3(signature = (path, *, rows=8192, bytes=4_194_304, options=None, context=None))]
pub(crate) fn open_structure_batches(
    py: Python<'_>,
    path: PathBuf,
    rows: usize,
    bytes: usize,
    options: Option<&crate::reading::PyReadOptions>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyStructureBatches> {
    let options = options.map_or_else(molframe::ReadOptions::new, |options| {
        options.native().clone()
    });
    let built = match context {
        Some(context) => context.build(context.token())?,
        None => ExecutionContext::default(),
    };
    let governing = &built;
    let source = py
        .detach(move || molframe::open_structure_batches(&path, &options, governing))
        .map_err(|error| crate::error::from_diagnostic(&molframe::Diagnostic::from(&error)))?;
    Ok(PyStructureBatches {
        source: Mutex::new(source),
        context: built,
        demand: BatchDemand::new(rows, bytes),
    })
}
