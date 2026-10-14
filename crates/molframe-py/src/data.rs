//! Manifest-driven table cache: list, verify and fetch, mechanically.
//!
//! The manifest, conversion, download and validation are `molframe::data`;
//! these functions pass paths in and rows out and release the GIL for the work.
//! `list` and `verify` read local files only; `fetch` is the one call that
//! touches the network, and only when the caller makes it. Errors keep their
//! registry codes (`E2201`, `E7101`, `E7102`, `E7901`).
//!
//! Cost: one manifest parse plus one cache-file read per entry; `fetch` adds
//! one download bounded by the caller's byte ceiling.

use molframe::chemistry::TableStore;
use molframe::data::{Fetched, Manifest, Verified};
use pyo3::prelude::*;
use std::path::PathBuf;

/// One manifest entry and the state of its cache file.
#[derive(Debug)]
#[pyclass(
    name = "DataEntry",
    frozen,
    skip_from_py_object,
    module = "molframe.data"
)]
struct PyDataEntry {
    inner: Verified,
}

#[pymethods]
impl PyDataEntry {
    #[getter]
    fn name(&self) -> &str {
        &self.inner.entry.name
    }

    #[getter]
    fn version(&self) -> &str {
        &self.inner.entry.version
    }

    #[getter]
    fn licence(&self) -> &str {
        &self.inner.entry.licence
    }

    #[getter]
    fn path(&self) -> PathBuf {
        self.inner.path.clone()
    }

    /// `missing`, `current`, `stale` or `invalid`.
    #[getter]
    fn status(&self) -> &'static str {
        self.inner.status.word()
    }

    /// Why the entry is stale or invalid; `None` otherwise.
    #[getter]
    fn detail(&self) -> Option<&str> {
        self.inner.status.detail()
    }

    fn __repr__(&self) -> String {
        format!(
            "DataEntry({}, {}, {})",
            self.inner.entry.name, self.inner.entry.version, self.inner.status
        )
    }
}

/// A table that `fetch` just stored.
#[derive(Debug)]
#[pyclass(
    name = "Fetched",
    frozen,
    skip_from_py_object,
    module = "molframe.data"
)]
struct PyFetched {
    inner: Fetched,
}

#[pymethods]
impl PyFetched {
    #[getter]
    fn path(&self) -> PathBuf {
        self.inner.path.clone()
    }

    #[getter]
    fn bytes(&self) -> usize {
        self.inner.bytes
    }

    #[getter]
    fn rows(&self) -> usize {
        self.inner.rows
    }

    fn __repr__(&self) -> String {
        format!(
            "Fetched({}, rows={})",
            self.inner.path.display(),
            self.inner.rows
        )
    }
}

fn entries(rows: Vec<Verified>) -> Vec<PyDataEntry> {
    rows.into_iter()
        .map(|inner| PyDataEntry { inner })
        .collect()
}

/// The directory a call works in: the caller's, else the store's own.
fn cache_directory(directory: Option<PathBuf>) -> PathBuf {
    match directory {
        Some(directory) => directory,
        None => TableStore::global().directory().to_path_buf(),
    }
}

/// Every entry of `manifest` with the state of its cache file; no network.
#[pyfunction]
#[pyo3(signature = (manifest, directory=None))]
fn list(
    py: Python<'_>,
    manifest: PathBuf,
    directory: Option<PathBuf>,
) -> PyResult<Vec<PyDataEntry>> {
    let rows = py.detach(move || -> PyResult<Vec<Verified>> {
        let manifest = Manifest::load(&manifest).map_err(crate::error::kernel)?;
        Ok(molframe::data::list(&manifest, &cache_directory(directory)))
    })?;
    Ok(entries(rows))
}

/// Re-reads the cache files of `manifest` (or the one `name`); no network.
#[pyfunction]
#[pyo3(signature = (manifest, directory=None, *, name=None))]
fn verify(
    py: Python<'_>,
    manifest: PathBuf,
    directory: Option<PathBuf>,
    name: Option<String>,
) -> PyResult<Vec<PyDataEntry>> {
    let rows = py.detach(move || -> PyResult<Vec<Verified>> {
        let manifest = Manifest::load(&manifest).map_err(crate::error::kernel)?;
        molframe::data::verify(&manifest, &cache_directory(directory), name.as_deref())
            .map_err(crate::error::kernel)
    })?;
    Ok(entries(rows))
}

/// The resource bounds of one download, all explicit.
#[derive(Clone, Copy, Debug)]
#[pyclass(
    name = "DownloadLimits",
    frozen,
    skip_from_py_object,
    module = "molframe.data"
)]
struct PyDownloadLimits {
    inner: molframe::adapters::DownloadOptions,
}

#[pymethods]
impl PyDownloadLimits {
    /// `max_bytes` bounds the body, `timeout_seconds` the whole request and
    /// `redirect_limit` the redirects followed (zero refuses them).
    #[new]
    #[pyo3(signature = (*, max_bytes, timeout_seconds, redirect_limit))]
    fn new(max_bytes: u64, timeout_seconds: u64, redirect_limit: usize) -> Self {
        Self {
            inner: molframe::adapters::DownloadOptions {
                max_bytes,
                timeout: std::time::Duration::from_secs(timeout_seconds),
                redirect_limit,
            },
        }
    }

    #[getter]
    fn max_bytes(&self) -> u64 {
        self.inner.max_bytes
    }

    #[getter]
    fn timeout_seconds(&self) -> u64 {
        self.inner.timeout.as_secs()
    }

    #[getter]
    fn redirect_limit(&self) -> usize {
        self.inner.redirect_limit
    }
}

/// Downloads `name`, converts and validates it, and stores it in the cache.
///
/// The only call that uses the network. The bytes are accepted only when their
/// SHA-256 matches the manifest; the resource bounds are the caller's.
#[pyfunction]
#[pyo3(signature = (manifest, name, limits, directory=None, *, replace=false))]
fn fetch(
    py: Python<'_>,
    manifest: PathBuf,
    name: String,
    limits: &PyDownloadLimits,
    directory: Option<PathBuf>,
    replace: bool,
) -> PyResult<PyFetched> {
    let options = limits.inner;
    let inner = py.detach(move || -> PyResult<Fetched> {
        let manifest = Manifest::load(&manifest).map_err(crate::error::kernel)?;
        molframe::data::fetch(
            &manifest,
            &name,
            &cache_directory(directory),
            replace,
            options,
        )
        .map_err(crate::error::kernel)
    })?;
    Ok(PyFetched { inner })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyDataEntry>()?;
    module.add_class::<PyFetched>()?;
    module.add_class::<PyDownloadLimits>()?;
    module.add_function(wrap_pyfunction!(list, module)?)?;
    module.add_function(wrap_pyfunction!(verify, module)?)?;
    module.add_function(wrap_pyfunction!(fetch, module)?)
}
