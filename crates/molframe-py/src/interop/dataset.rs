//! Manifest-backed datasets: filter, split and batch without loading coordinates.

use crate::bindings::{PyStructure, findings_error};
use molframe::interop::{
    Dataset, DatasetFilter, DatasetSplit, LoadError, ManifestEntry, SplitOptions, SplitRatios,
    SplitStrategy,
};
use pyo3::prelude::*;
use std::collections::BTreeMap;
use std::path::PathBuf;

/// One structure of a dataset and what is known of it without loading it.
#[derive(Clone, Debug)]
#[pyclass(
    name = "ManifestEntry",
    frozen,
    skip_from_py_object,
    module = "molframe.interop"
)]
pub(crate) struct PyManifestEntry {
    inner: ManifestEntry,
}

#[pymethods]
impl PyManifestEntry {
    #[new]
    #[pyo3(signature = (
        id,
        path,
        atom_count,
        *,
        resolution=None,
        method=None,
        deposition_date=None,
        sequence=None,
        structure_cluster=None,
        tags=Vec::new(),
        statistics=BTreeMap::new(),
    ))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        id: &str,
        path: PathBuf,
        atom_count: u64,
        resolution: Option<f32>,
        method: Option<&str>,
        deposition_date: Option<&str>,
        sequence: Option<&str>,
        structure_cluster: Option<&str>,
        tags: Vec<String>,
        statistics: BTreeMap<String, f64>,
    ) -> Self {
        Self {
            inner: ManifestEntry {
                id: id.into(),
                path,
                atom_count,
                resolution,
                method: method.map(Into::into),
                deposition_date: deposition_date.map(Into::into),
                sequence: sequence.map(Into::into),
                structure_cluster: structure_cluster.map(Into::into),
                tags: tags.into_iter().map(Into::into).collect(),
                statistics: statistics
                    .into_iter()
                    .map(|(name, value)| (name.into(), value))
                    .collect(),
            },
        }
    }

    #[getter]
    fn id(&self) -> String {
        self.inner.id.to_string()
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
    fn method(&self) -> Option<String> {
        self.inner.method.as_deref().map(str::to_owned)
    }

    #[getter]
    fn deposition_date(&self) -> Option<String> {
        self.inner.deposition_date.as_deref().map(str::to_owned)
    }

    #[getter]
    fn sequence(&self) -> Option<String> {
        self.inner.sequence.as_deref().map(str::to_owned)
    }

    #[getter]
    fn structure_cluster(&self) -> Option<String> {
        self.inner.structure_cluster.as_deref().map(str::to_owned)
    }

    #[getter]
    fn tags(&self) -> Vec<String> {
        self.inner.tags.iter().map(ToString::to_string).collect()
    }

    #[getter]
    fn statistics(&self) -> BTreeMap<String, f64> {
        self.inner
            .statistics
            .iter()
            .map(|(name, value)| (name.to_string(), *value))
            .collect()
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }

    fn __repr__(&self) -> String {
        format!(
            "ManifestEntry(id={:?}, atoms={})",
            self.inner.id, self.inner.atom_count
        )
    }
}

/// The three partitions of a split, and what the policy warns about.
#[derive(Debug)]
#[pyclass(
    name = "DatasetSplit",
    frozen,
    skip_from_py_object,
    module = "molframe.interop"
)]
pub(crate) struct PyDatasetSplit {
    inner: DatasetSplit,
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

    /// What the strategy cannot guarantee, in words.
    #[getter]
    fn warnings(&self) -> Vec<String> {
        self.inner
            .warnings
            .iter()
            .map(ToString::to_string)
            .collect()
    }

    fn __repr__(&self) -> String {
        format!(
            "DatasetSplit(train={}, validation={}, test={})",
            self.inner.train.len(),
            self.inner.validation.len(),
            self.inner.test.len()
        )
    }
}

/// An immutable, lazy list of structures described by a manifest.
///
/// Filtering, splitting and batching read only the manifest; a structure is
/// loaded when `load` asks for it.
#[derive(Clone, Debug)]
#[pyclass(
    name = "Dataset",
    frozen,
    skip_from_py_object,
    module = "molframe.interop"
)]
pub(crate) struct PyDataset {
    inner: Dataset,
}

fn at(dataset: &Dataset, index: isize) -> PyResult<usize> {
    let length = dataset.len();
    let resolved = if index < 0 {
        length.checked_sub(index.unsigned_abs())
    } else {
        Some(index.cast_unsigned())
    };
    match resolved {
        Some(position) if position < length => Ok(position),
        _ => Err(crate::error::index("dataset index is out of range")),
    }
}

#[pymethods]
impl PyDataset {
    /// A dataset of the given entries; identifiers must be unique.
    #[new]
    fn new(entries: Vec<PyRef<'_, PyManifestEntry>>) -> PyResult<Self> {
        let entries = entries
            .into_iter()
            .map(|entry| entry.inner.clone())
            .collect();
        Dataset::new(entries)
            .map(|inner| Self { inner })
            .map_err(crate::error::kernel)
    }

    /// Reads a JSON manifest: a list of entries, or an object with `entries`.
    #[staticmethod]
    fn from_manifest(py: Python<'_>, path: PathBuf) -> PyResult<Self> {
        py.detach(move || Dataset::from_manifest(&path))
            .map(|inner| Self { inner })
            .map_err(crate::error::kernel)
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __getitem__(&self, index: isize) -> PyResult<PyManifestEntry> {
        let position = at(&self.inner, index)?;
        self.inner
            .entries()
            .nth(position)
            .map(|entry| PyManifestEntry {
                inner: entry.clone(),
            })
            .ok_or_else(|| crate::error::index("dataset index is out of range"))
    }

    /// The selected entries in manifest order.
    #[getter]
    fn entries(&self) -> Vec<PyManifestEntry> {
        self.inner
            .entries()
            .map(|entry| PyManifestEntry {
                inner: entry.clone(),
            })
            .collect()
    }

    /// The entries that satisfy every stated criterion.
    #[pyo3(signature = (
        *,
        resolution_below=None,
        method=None,
        minimum_atoms=None,
        maximum_atoms=None,
        tag=None,
    ))]
    fn filter(
        &self,
        resolution_below: Option<f32>,
        method: Option<&str>,
        minimum_atoms: Option<u64>,
        maximum_atoms: Option<u64>,
        tag: Option<&str>,
    ) -> PyResult<Self> {
        let filter = DatasetFilter {
            resolution_below,
            method: method.map(Into::into),
            minimum_atoms,
            maximum_atoms,
            tag: tag.map(Into::into),
        };
        self.inner
            .filter(&filter)
            .map(|inner| Self { inner })
            .map_err(crate::error::kernel)
    }

    /// Train, validation and test partitions under one stated strategy.
    ///
    /// `strategy` is `"sequence_identity"` (with `threshold`), `"random"`
    /// (with `seed`), `"structural_cluster"` or `"temporal"`. The ratios are
    /// stated, not defaulted, and must sum to one.
    #[pyo3(signature = (strategy, *, ratios, threshold=None, seed=None))]
    fn split(
        &self,
        py: Python<'_>,
        strategy: &str,
        ratios: (f64, f64, f64),
        threshold: Option<f64>,
        seed: Option<u64>,
    ) -> PyResult<PyDatasetSplit> {
        let options = SplitOptions {
            strategy: SplitStrategy::from_parts(strategy, threshold, seed)
                .map_err(crate::error::kernel)?,
            ratios: SplitRatios::new(ratios.0, ratios.1, ratios.2).map_err(crate::error::kernel)?,
        };
        let dataset = self.inner.clone();
        py.detach(move || dataset.split(&options))
            .map(|inner| PyDatasetSplit { inner })
            .map_err(crate::error::kernel)
    }

    /// Consecutive batches of at most `size` entries.
    fn batches(&self, size: usize) -> PyResult<Vec<Self>> {
        self.inner
            .batches(size)
            .map(|batches| batches.into_iter().map(|inner| Self { inner }).collect())
            .map_err(crate::error::kernel)
    }

    /// Reads the structure of entry `index`.
    fn load(&self, py: Python<'_>, index: isize) -> PyResult<PyStructure> {
        let position = at(&self.inner, index)?;
        let dataset = self.inner.clone();
        py.detach(move || dataset.load_with(position, |entry| molframe::read(&entry.path)))
            .map(PyStructure::new)
            .map_err(|error| match error {
                LoadError::Dataset(error) => crate::error::kernel(error),
                LoadError::Loader(findings) => findings_error(&findings),
            })
    }

    fn __repr__(&self) -> String {
        format!("Dataset(entries={})", self.inner.len())
    }
}
