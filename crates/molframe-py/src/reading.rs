//! Reading structures: the options a read takes, and the readers.

use crate::bindings::{PyStructure, findings_error};
use crate::execution::PyExecutionContext;
use molframe::{
    AmbiguousResidueBoundaryPolicy, CategoryFilter, Diagnostic, Format, Limits,
    MissingElementPolicy, ParseMode, ReadOptions,
};
use pyo3::prelude::*;
use pyo3::pybacked::PyBackedBytes;
use pyo3::types::PyList;
use std::fmt;
use std::path::PathBuf;

/// Parses one word of the read vocabulary through Rust.
fn parse<T>(value: &str) -> PyResult<T>
where
    T: std::str::FromStr<Err = molframe::PolicyParseError>,
{
    value.parse().map_err(crate::error::kernel)
}

/// `value`, or `default` when none was stated.
fn or_default<T>(value: Option<T>, default: T) -> T {
    match value {
        Some(value) => value,
        None => default,
    }
}

/// Everything a read decides, stated once and reusable across reads.
///
/// Omitted arguments keep the library's choices: automatic format detection,
/// permissive parsing, every model, every category, and the default input
/// limits. `only_categories` and `skip_categories` choose which optional
/// mmCIF/BinaryCIF categories or PDB records are read; what a structure cannot
/// be built without is always read.
#[derive(Clone, Debug)]
#[pyclass(name = "ReadOptions", frozen, skip_from_py_object, module = "molframe")]
pub(crate) struct PyReadOptions {
    inner: ReadOptions,
}

impl PyReadOptions {
    /// The native options this value stands for.
    pub(crate) const fn native(&self) -> &ReadOptions {
        &self.inner
    }
}

#[pymethods]
impl PyReadOptions {
    #[new]
    #[pyo3(signature = (
        *,
        format="auto",
        mode="permissive",
        first_model_only=false,
        coordinates_only=false,
        discard_hydrogens=false,
        missing_element="preserve_unknown",
        ambiguous_residue_boundary="reject",
        only_categories=None,
        skip_categories=None,
        max_decompressed_bytes=None,
        max_compression_ratio=None,
        max_rows_per_category=None,
        max_nesting_depth=None,
        max_dictionary_entries=None,
    ))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        format: &str,
        mode: &str,
        first_model_only: bool,
        coordinates_only: bool,
        discard_hydrogens: bool,
        missing_element: &str,
        ambiguous_residue_boundary: &str,
        only_categories: Option<Vec<String>>,
        skip_categories: Option<Vec<String>>,
        max_decompressed_bytes: Option<u64>,
        max_compression_ratio: Option<u64>,
        max_rows_per_category: Option<u64>,
        max_nesting_depth: Option<u32>,
        max_dictionary_entries: Option<u32>,
    ) -> PyResult<Self> {
        let categories = match (only_categories, skip_categories) {
            (None, None) => CategoryFilter::All,
            (Some(names), None) => CategoryFilter::only(names),
            (None, Some(names)) => CategoryFilter::except(names),
            (Some(_), Some(_)) => {
                return Err(crate::error::value(
                    "pass only_categories or skip_categories, not both",
                ));
            }
        };
        let defaults = Limits::default();
        let limits = Limits {
            decompressed_bytes: or_default(max_decompressed_bytes, defaults.decompressed_bytes),
            compression_ratio: or_default(max_compression_ratio, defaults.compression_ratio),
            rows_per_category: or_default(max_rows_per_category, defaults.rows_per_category),
            nesting_depth: or_default(max_nesting_depth, defaults.nesting_depth),
            dictionary_entries: or_default(max_dictionary_entries, defaults.dictionary_entries),
        };
        Ok(Self {
            inner: ReadOptions::new()
                .format(parse::<Format>(format)?)
                .mode(parse::<ParseMode>(mode)?)
                .only_first_model(first_model_only)
                .only_atomic_coords(coordinates_only)
                .discard_hydrogens(discard_hydrogens)
                .missing_element_policy(parse::<MissingElementPolicy>(missing_element)?)
                .ambiguous_residue_boundary_policy(parse::<AmbiguousResidueBoundaryPolicy>(
                    ambiguous_residue_boundary,
                )?)
                .categories(categories)
                .limits(limits),
        })
    }

    #[getter]
    fn format(&self) -> String {
        crate::policy::snake(self.inner.format.name())
    }

    #[getter]
    fn mode(&self) -> String {
        crate::policy::snake(self.inner.mode.name())
    }

    #[getter]
    const fn first_model_only(&self) -> bool {
        self.inner.only_first_model
    }

    #[getter]
    const fn coordinates_only(&self) -> bool {
        self.inner.only_atomic_coords
    }

    #[getter]
    const fn discard_hydrogens(&self) -> bool {
        self.inner.discard_hydrogens
    }

    #[getter]
    fn missing_element(&self) -> String {
        crate::policy::snake(self.inner.missing_element_policy.name())
    }

    #[getter]
    fn ambiguous_residue_boundary(&self) -> String {
        crate::policy::snake(self.inner.ambiguous_residue_boundary_policy.name())
    }

    #[getter]
    fn only_categories(&self) -> Option<Vec<String>> {
        match &self.inner.categories {
            CategoryFilter::Only(names) => Some(names.iter().map(ToString::to_string).collect()),
            _ => None,
        }
    }

    #[getter]
    fn skip_categories(&self) -> Option<Vec<String>> {
        match &self.inner.categories {
            CategoryFilter::Except(names) => Some(names.iter().map(ToString::to_string).collect()),
            _ => None,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "ReadOptions(format={}, mode={})",
            self.format(),
            self.mode()
        )
    }
}

/// The options a call reads under: `options`, with `format` standing in for the
/// format it leaves to detection.
fn resolve(options: Option<&PyReadOptions>, format: Option<&str>) -> PyResult<ReadOptions> {
    let mut resolved = options.map_or_else(ReadOptions::new, |options| options.inner.clone());
    if let Some(name) = format {
        if resolved.format != Format::Auto {
            return Err(crate::error::value(
                "the format is stated twice: pass format= or ReadOptions(format=), not both",
            ));
        }
        resolved = resolved.format(parse::<Format>(name)?);
    }
    Ok(resolved)
}

struct PythonBytes(PyBackedBytes);

impl AsRef<[u8]> for PythonBytes {
    fn as_ref(&self) -> &[u8] {
        self.0.as_ref()
    }
}

impl fmt::Debug for PythonBytes {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("PythonBytes")
            .field("len", &self.0.len())
            .finish()
    }
}

/// What a read produces: the structure and what was wrong with its file.
type Loaded = (PyStructure, Vec<Diagnostic>);

fn load_buffer(
    py: Python<'_>,
    input: &molframe::InputBuffer,
    name: Option<&str>,
    options: &ReadOptions,
    context: Option<&PyExecutionContext>,
) -> PyResult<Loaded> {
    let input = input.clone();
    let name = name.map(str::to_owned);
    let options = options.clone();
    crate::execution::run(py, context, move |context| {
        molframe::read_buffer_in(&input, name.as_deref(), &options, context)
    })?
    .map(|(structure, findings)| (PyStructure::new(structure), findings))
    .map_err(|findings| findings_error(&findings))
}

fn load(
    py: Python<'_>,
    source: &Bound<'_, PyAny>,
    name: Option<&str>,
    options: &ReadOptions,
    context: Option<&PyExecutionContext>,
) -> PyResult<Loaded> {
    if let Ok(path) = source.extract::<PathBuf>() {
        let options = options.clone();
        return crate::execution::run(py, context, move |context| {
            molframe::read_with_options_in(&path, &options, context)
        })?
        .map(|(structure, findings)| (PyStructure::new(structure), findings))
        .map_err(|findings| findings_error(&findings));
    }
    let bytes = source.extract::<PyBackedBytes>()?;
    let input = molframe::InputBuffer::from_owner(PythonBytes(bytes));
    load_buffer(py, &input, name, options, context)
}

/// Named in-memory bytes, read without copying them.
#[derive(Clone, Debug)]
#[pyclass(name = "Reader", frozen, skip_from_py_object)]
pub(crate) struct PyReader {
    input: molframe::InputBuffer,
    name: Option<Box<str>>,
}

#[pymethods]
impl PyReader {
    #[new]
    #[pyo3(signature = (data, *, name=None))]
    fn new(data: PyBackedBytes, name: Option<&str>) -> Self {
        Self {
            input: molframe::InputBuffer::from_owner(PythonBytes(data)),
            name: name.map(Into::into),
        }
    }

    #[pyo3(signature = (*, format=None, options=None, context=None))]
    fn read(
        &self,
        py: Python<'_>,
        format: Option<&str>,
        options: Option<&PyReadOptions>,
        context: Option<&PyExecutionContext>,
    ) -> PyResult<PyStructure> {
        let options = resolve(options, format)?;
        load_buffer(py, &self.input, self.name.as_deref(), &options, context)
            .map(|(structure, _)| structure)
    }

    #[getter]
    fn byte_length(&self) -> usize {
        self.input.len()
    }
}

/// Reads a path or bytes into a structure, discarding what was wrong with the file.
///
/// Use `read_with_diagnostics` to keep the findings. `format` names the format
/// when neither the content nor the file name settles it; `options` states
/// every other decision; `context` bounds the work and lets it be cancelled.
#[pyfunction]
#[pyo3(signature = (source, *, name=None, format=None, options=None, context=None))]
pub(crate) fn read(
    py: Python<'_>,
    source: &Bound<'_, PyAny>,
    name: Option<&str>,
    format: Option<&str>,
    options: Option<&PyReadOptions>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyStructure> {
    let options = resolve(options, format)?;
    load(py, source, name, &options, context).map(|(structure, _)| structure)
}

/// Reads a structure and returns, with it, every finding about the file.
///
/// A file that parses is not the same as a file that is right: the findings are
/// `Diagnostic` values with a code, a message, a remedy and the byte span they
/// concern, ordered as the reader found them.
#[pyfunction]
#[pyo3(signature = (source, *, name=None, format=None, options=None, context=None))]
pub(crate) fn read_with_diagnostics<'py>(
    py: Python<'py>,
    source: &Bound<'py, PyAny>,
    name: Option<&str>,
    format: Option<&str>,
    options: Option<&PyReadOptions>,
    context: Option<&PyExecutionContext>,
) -> PyResult<(PyStructure, Bound<'py, PyList>)> {
    let options = resolve(options, format)?;
    let (structure, findings) = load(py, source, name, &options, context)?;
    Ok((structure, crate::error::diagnostic_list(py, &findings)?))
}
