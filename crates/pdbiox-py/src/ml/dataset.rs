//! Declarative Python views over the native lazy dataset engine.

use crate::errors::read_error;
use crate::structure::PyStructure;
use pdbiox::{
    Dataset, DatasetFilter, DatasetSplit, DatasetWarning, ManifestEntry, SplitOptions, SplitRatios,
    SplitStrategy,
};
use pyo3::exceptions::{PyIndexError, PyRuntimeError};
use pyo3::prelude::*;
use pyo3::types::PyList;
use std::collections::BTreeMap;
use std::path::PathBuf;

#[pyclass(name = "DatasetEntry", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyDatasetEntry {
    inner: ManifestEntry,
}

#[pymethods]
impl PyDatasetEntry {
    #[new]
    #[pyo3(signature = (id, path, atom_count, *, resolution=None, method=None, deposition_date=None, sequence=None, structure_cluster=None, tags=Vec::new(), statistics=None))]
    fn new(
        id: String,
        path: PathBuf,
        atom_count: u64,
        resolution: Option<f32>,
        method: Option<String>,
        deposition_date: Option<String>,
        sequence: Option<String>,
        structure_cluster: Option<String>,
        tags: Vec<String>,
        statistics: Option<BTreeMap<String, f64>>,
    ) -> Self {
        Self {
            inner: ManifestEntry {
                id: id.into_boxed_str(),
                path,
                atom_count,
                resolution,
                method: method.map(String::into_boxed_str),
                deposition_date: deposition_date.map(String::into_boxed_str),
                sequence: sequence.map(String::into_boxed_str),
                structure_cluster: structure_cluster.map(String::into_boxed_str),
                tags: tags.into_iter().map(String::into_boxed_str).collect(),
                statistics: match statistics {
                    Some(values) => values
                        .into_iter()
                        .map(|(name, value)| (name.into_boxed_str(), value))
                        .collect(),
                    None => BTreeMap::new(),
                },
            },
        }
    }

    #[getter]
    fn id(&self) -> &str {
        &self.inner.id
    }

    #[getter]
    fn path(&self) -> PathBuf {
        self.inner.path.clone()
    }

    #[getter]
    const fn atom_count(&self) -> u64 {
        self.inner.atom_count
    }

    #[getter]
    const fn resolution(&self) -> Option<f32> {
        self.inner.resolution
    }

    #[getter]
    fn method(&self) -> Option<&str> {
        self.inner.method.as_deref()
    }

    #[getter]
    fn deposition_date(&self) -> Option<&str> {
        self.inner.deposition_date.as_deref()
    }

    #[getter]
    fn sequence(&self) -> Option<&str> {
        self.inner.sequence.as_deref()
    }

    #[getter]
    fn structure_cluster(&self) -> Option<&str> {
        self.inner.structure_cluster.as_deref()
    }

    #[getter]
    fn tags(&self) -> Vec<&str> {
        self.inner.tags.iter().map(AsRef::as_ref).collect()
    }

    #[getter]
    fn statistics(&self) -> BTreeMap<&str, f64> {
        self.inner
            .statistics
            .iter()
            .map(|(name, value)| (name.as_ref(), *value))
            .collect()
    }

    fn __repr__(&self) -> String {
        format!(
            "DatasetEntry(id={:?}, path={}, atom_count={})",
            self.inner.id,
            self.inner.path.display(),
            self.inner.atom_count
        )
    }
}

#[pyclass(name = "DatasetFilter", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyDatasetFilter {
    pub(super) inner: DatasetFilter,
}

#[pymethods]
impl PyDatasetFilter {
    #[new]
    #[pyo3(signature = (*, resolution_below=None, method=None, minimum_atoms=None, maximum_atoms=None, tag=None))]
    fn new(
        resolution_below: Option<f32>,
        method: Option<String>,
        minimum_atoms: Option<u64>,
        maximum_atoms: Option<u64>,
        tag: Option<String>,
    ) -> Self {
        Self {
            inner: DatasetFilter {
                resolution_below,
                method: method.map(Into::into),
                minimum_atoms,
                maximum_atoms,
                tag: tag.map(Into::into),
            },
        }
    }

    #[getter]
    const fn resolution_below(&self) -> Option<f32> {
        self.inner.resolution_below
    }

    #[getter]
    fn method(&self) -> Option<&str> {
        self.inner.method.as_deref()
    }

    #[getter]
    const fn minimum_atoms(&self) -> Option<u64> {
        self.inner.minimum_atoms
    }

    #[getter]
    const fn maximum_atoms(&self) -> Option<u64> {
        self.inner.maximum_atoms
    }

    #[getter]
    fn tag(&self) -> Option<&str> {
        self.inner.tag.as_deref()
    }
}

#[pyclass(name = "SplitRatios", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PySplitRatios {
    pub(super) inner: SplitRatios,
}

#[pymethods]
impl PySplitRatios {
    #[new]
    fn new(train: f64, validation: f64, test: f64) -> PyResult<Self> {
        SplitRatios::new(train, validation, test)
            .map(|inner| Self { inner })
            .map_err(crate::errors::dataset_error)
    }

    #[getter]
    const fn train(&self) -> f64 {
        self.inner.train
    }

    #[getter]
    const fn validation(&self) -> f64 {
        self.inner.validation
    }

    #[getter]
    const fn test(&self) -> f64 {
        self.inner.test
    }

    fn __repr__(&self) -> String {
        format!(
            "SplitRatios(train={}, validation={}, test={})",
            self.inner.train, self.inner.validation, self.inner.test
        )
    }
}

#[pyclass(name = "SplitOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PySplitOptions {
    pub(super) inner: SplitOptions,
}

#[pymethods]
impl PySplitOptions {
    #[new]
    #[pyo3(signature = (strategy, *, ratios=None, train=0.8, validation=0.1, test=0.1))]
    fn new(
        strategy: PySplitStrategy,
        ratios: Option<PySplitRatios>,
        train: f64,
        validation: f64,
        test: f64,
    ) -> PyResult<Self> {
        let ratios = ratios
            .map_or_else(
                || SplitRatios::new(train, validation, test),
                |ratios| Ok(ratios.inner),
            )
            .map_err(crate::errors::dataset_error)?;
        Ok(Self {
            inner: SplitOptions {
                strategy: strategy.inner,
                ratios,
            },
        })
    }

    #[getter]
    fn strategy(&self) -> PySplitStrategy {
        PySplitStrategy {
            inner: self.inner.strategy,
        }
    }

    #[getter]
    const fn ratios(&self) -> PySplitRatios {
        PySplitRatios {
            inner: self.inner.ratios,
        }
    }
}

#[pyclass(name = "SplitStrategy", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PySplitStrategy {
    inner: SplitStrategy,
}

#[pymethods]
impl PySplitStrategy {
    #[staticmethod]
    fn sequence_identity(threshold: f64) -> Self {
        Self {
            inner: SplitStrategy::SequenceIdentity { threshold },
        }
    }

    #[staticmethod]
    fn structural_cluster() -> Self {
        Self {
            inner: SplitStrategy::StructuralCluster,
        }
    }

    #[staticmethod]
    fn temporal() -> Self {
        Self {
            inner: SplitStrategy::Temporal,
        }
    }

    #[staticmethod]
    fn random(seed: u64) -> Self {
        Self {
            inner: SplitStrategy::Random { seed },
        }
    }
}

#[pyclass(name = "DatasetWarning", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyDatasetWarning {
    RandomSplitMayLeak,
}

#[pyclass(name = "Dataset", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyDataset {
    inner: Dataset,
}

#[pymethods]
impl PyDataset {
    #[new]
    fn new(path: PathBuf) -> PyResult<Self> {
        let result = Dataset::from_manifest(&path)
            .map(|inner| Self { inner })
            .map_err(dataset_error);
        drop(path);
        result
    }

    #[staticmethod]
    fn from_entries(entries: &Bound<'_, PyList>) -> PyResult<Self> {
        let entries = entries
            .iter()
            .map(|entry| {
                entry
                    .extract::<PyRef<'_, PyDatasetEntry>>()
                    .map_err(pyo3::PyErr::from)
                    .map(|entry| entry.inner.clone())
            })
            .collect::<PyResult<Vec<_>>>()?;
        Dataset::new(entries)
            .map(|inner| Self { inner })
            .map_err(dataset_error)
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    fn __getitem__(&self, index: isize) -> PyResult<PyDatasetEntry> {
        let index = normalise_index(index, self.inner.len())?;
        self.inner
            .entries()
            .nth(index)
            .cloned()
            .map(|inner| PyDatasetEntry { inner })
            .ok_or_else(|| PyIndexError::new_err(index))
    }

    fn load(&self, py: Python<'_>, index: isize) -> PyResult<PyStructure> {
        let index = normalise_index(index, self.inner.len())?;
        match self
            .inner
            .load_with(index, |entry| pdbiox::read(&entry.path))
        {
            Ok(structure) => Ok(PyStructure::new(structure)),
            Err(pdbiox::LoadError::Dataset(error)) => Err(dataset_error(error)),
            Err(pdbiox::LoadError::Loader(findings)) => Err(read_error(py, &findings)),
        }
    }

    #[pyo3(signature = (*, resolution_below=None, method=None, minimum_atoms=None, maximum_atoms=None, tag=None))]
    fn filter(
        &self,
        resolution_below: Option<f32>,
        method: Option<String>,
        minimum_atoms: Option<u64>,
        maximum_atoms: Option<u64>,
        tag: Option<String>,
    ) -> PyResult<Self> {
        self.filter_spec(&PyDatasetFilter::new(
            resolution_below,
            method,
            minimum_atoms,
            maximum_atoms,
            tag,
        ))
    }

    fn filter_spec(&self, filter: &PyDatasetFilter) -> PyResult<Self> {
        self.inner
            .filter(&filter.inner)
            .map(|inner| Self { inner })
            .map_err(dataset_error)
    }

    #[pyo3(signature = (strategy, *, train=0.8, validation=0.1, test=0.1))]
    fn split(
        &self,
        strategy: PySplitStrategy,
        train: f64,
        validation: f64,
        test: f64,
    ) -> PyResult<PyDatasetSplit> {
        self.split_options(PySplitOptions::new(
            strategy, None, train, validation, test,
        )?)
    }

    fn split_options(&self, options: PySplitOptions) -> PyResult<PyDatasetSplit> {
        self.inner
            .split(&options.inner)
            .map(PyDatasetSplit::from)
            .map_err(dataset_error)
    }

    fn batches(&self, batch_size: usize) -> PyResult<Vec<Self>> {
        self.inner
            .batches(batch_size)
            .map(|batches| batches.into_iter().map(|inner| Self { inner }).collect())
            .map_err(dataset_error)
    }
}

#[pyclass(name = "DatasetSplit", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyDatasetSplit {
    inner: DatasetSplit,
}

impl From<DatasetSplit> for PyDatasetSplit {
    fn from(inner: DatasetSplit) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PyDatasetSplit {
    #[getter]
    fn train(&self) -> PyDataset {
        PyDataset {
            inner: self.inner.train.clone(),
        }
    }

    #[getter]
    fn validation(&self) -> PyDataset {
        PyDataset {
            inner: self.inner.validation.clone(),
        }
    }

    #[getter]
    fn test(&self) -> PyDataset {
        PyDataset {
            inner: self.inner.test.clone(),
        }
    }

    #[getter]
    fn warnings(&self) -> PyResult<Vec<PyDatasetWarning>> {
        self.inner
            .warnings
            .iter()
            .copied()
            .map(|warning| match warning {
                DatasetWarning::RandomSplitMayLeak => Ok(PyDatasetWarning::RandomSplitMayLeak),
                _ => Err(PyRuntimeError::new_err(
                    "dataset warning is not exposed by this binding version",
                )),
            })
            .collect()
    }
}

fn normalise_index(index: isize, length: usize) -> PyResult<usize> {
    let length = isize::try_from(length).map_err(|_| PyIndexError::new_err(index))?;
    let index = if index < 0 { length + index } else { index };
    if index < 0 || index >= length {
        Err(PyIndexError::new_err(index))
    } else {
        usize::try_from(index).map_err(|_| PyIndexError::new_err(index))
    }
}

fn dataset_error(error: pdbiox::DatasetError) -> PyErr {
    crate::errors::dataset_error(error)
}
