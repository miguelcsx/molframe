//! Native k-mer index and spaced-seed adapters.

use ::pdbiox;
use pyo3::create_exception;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

create_exception!(_native, KmerTableError, PyValueError);
create_exception!(_native, SeedPatternError, PyValueError);

#[pyclass(name = "KmerStorage", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyKmerStorage(pub(crate) pdbiox::seq::KmerStorage);

#[pymethods]
impl PyKmerStorage {
    #[staticmethod]
    fn exact() -> Self {
        Self(pdbiox::seq::KmerStorage::Exact)
    }

    #[staticmethod]
    fn bucketed(buckets: usize) -> PyResult<Self> {
        if buckets == 0 {
            return Err(KmerTableError::new_err(
                "bucketed k-mer table needs buckets",
            ));
        }
        Ok(Self(pdbiox::seq::KmerStorage::Bucketed { buckets }))
    }

    #[getter]
    fn buckets(&self) -> Option<usize> {
        match self.0 {
            pdbiox::seq::KmerStorage::Exact => None,
            pdbiox::seq::KmerStorage::Bucketed { buckets } => Some(buckets),
        }
    }

    #[getter]
    fn is_exact(&self) -> bool {
        matches!(self.0, pdbiox::seq::KmerStorage::Exact)
    }
}

#[pyclass(name = "SeedPattern", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySeedPattern(pub(crate) pdbiox::seq::SeedPattern);

#[pymethods]
impl PySeedPattern {
    #[new]
    fn new(mask: Vec<bool>) -> PyResult<Self> {
        pdbiox::seq::SeedPattern::new(&mask)
            .map(Self)
            .map_err(|error| SeedPatternError::new_err(error.to_string()))
    }

    #[getter]
    fn span(&self) -> usize {
        self.0.span()
    }

    #[getter]
    fn weight(&self) -> usize {
        self.0.weight()
    }
}

#[pyclass(name = "KmerTableOptions", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyKmerTableOptions(pub(crate) pdbiox::seq::KmerTableOptions);

#[pymethods]
impl PyKmerTableOptions {
    #[new]
    #[pyo3(signature = (k, *, pattern=None, storage=None))]
    fn new(k: usize, pattern: Option<PySeedPattern>, storage: Option<PyKmerStorage>) -> Self {
        Self(pdbiox::seq::KmerTableOptions {
            k,
            pattern: pattern.map(|value| value.0),
            storage: storage.map_or(pdbiox::seq::KmerStorage::Exact, |value| value.0),
        })
    }

    #[getter]
    fn k(&self) -> usize {
        self.0.k
    }

    #[getter]
    fn pattern(&self) -> Option<PySeedPattern> {
        self.0.pattern.clone().map(PySeedPattern)
    }

    #[getter]
    fn storage(&self) -> PyKmerStorage {
        PyKmerStorage(self.0.storage)
    }
}

#[pyclass(name = "KmerHit", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyKmerHit {
    #[pyo3(get)]
    sequence: usize,
    #[pyo3(get)]
    position: usize,
}

#[pyclass(name = "KmerTable", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyKmerTable(pub(crate) pdbiox::seq::KmerTable);

#[pymethods]
impl PyKmerTable {
    #[staticmethod]
    fn build(sequences: Vec<Vec<u8>>, options: &PyKmerTableOptions) -> PyResult<Self> {
        let borrowed = sequences.iter().map(Vec::as_slice).collect::<Vec<_>>();
        pdbiox::seq::KmerTable::build(&borrowed, &options.0)
            .map(Self)
            .map_err(|error| KmerTableError::new_err(error.to_string()))
    }

    fn query(&self, word: Vec<u8>) -> Vec<PyKmerHit> {
        self.0
            .query(&word)
            .iter()
            .copied()
            .map(|hit| PyKmerHit {
                sequence: hit.sequence,
                position: hit.position,
            })
            .collect()
    }

    #[getter]
    fn word_length(&self) -> usize {
        self.0.word_length()
    }
}

impl From<pdbiox::seq::KmerHit> for PyKmerHit {
    fn from(value: pdbiox::seq::KmerHit) -> Self {
        Self {
            sequence: value.sequence,
            position: value.position,
        }
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("KmerTableError", module.py().get_type::<KmerTableError>())?;
    module.add(
        "SeedPatternError",
        module.py().get_type::<SeedPatternError>(),
    )?;
    module.add_class::<PyKmerStorage>()?;
    module.add_class::<PySeedPattern>()?;
    module.add_class::<PyKmerTableOptions>()?;
    module.add_class::<PyKmerHit>()?;
    module.add_class::<PyKmerTable>()?;
    Ok(())
}
