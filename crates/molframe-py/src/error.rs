//! Python exceptions for `MolFrame` failures.
//!
//! The classes live in `molframe.errors` — plain Python, so each can inherit
//! the built-in it should be catchable as — and this module only chooses which
//! one a diagnostic becomes and fills its structured fields. Nothing here
//! decides what went wrong; the kernels name a registered code and the class
//! follows from the code.

use molframe::{Class, Code, Diagnostic};
use pyo3::prelude::*;
use pyo3::sync::PyOnceLock;
use pyo3::types::{PyDict, PyTuple};

/// The exception class a diagnostic code is raised as.
///
/// The exact codes come first because they name a failure Python code wants to
/// tell apart; everything else follows its class, the thousands digit of the
/// code.
pub(crate) const fn class_name(code: Code) -> &'static str {
    // Compared by number: `Code` is a plain value and these are its few
    // exceptions to the by-class rule.
    match code.number() {
        1902 | 7001 => "MemoryBudgetError",
        1904 => "Cancelled",
        7101 | 7901 => "MolframeIOError",
        1901 | 1903 => "ResourceError",
        _ => match code.class() {
            Class::Syntax => "ParseError",
            Class::Schema => "SchemaError",
            Class::Consistency => "ConsistencyError",
            Class::Conversion => "ConversionError",
            Class::Geometry => "GeometryError",
            Class::Policy => "PolicyError",
            Class::Resource => "ResourceError",
            // Internal, and any class a newer library adds that this table does
            // not know: neither is something a caller can act on.
            Class::Internal | _ => "InternalError",
        },
    }
}

/// Every exception class name this module can raise, for the registration test.
#[cfg(test)]
pub(crate) const CLASS_NAMES: &[&str] = &[
    "MolframeError",
    "MolframeValueError",
    "MolframeKeyError",
    "MolframeIndexError",
    "MolframeTypeError",
    "ParseError",
    "SchemaError",
    "ConsistencyError",
    "ConversionError",
    "QueryError",
    "GeometryError",
    "PolicyError",
    "ResourceError",
    "MolframeIOError",
    "MemoryBudgetError",
    "Cancelled",
    "InternalError",
    "MolframeWarning",
    "QueryWarning",
    "Diagnostic",
];

static ERRORS: PyOnceLock<Py<PyModule>> = PyOnceLock::new();

fn errors_module(py: Python<'_>) -> PyResult<&Bound<'_, PyModule>> {
    Ok(ERRORS
        .get_or_try_init(py, || py.import("molframe.errors").map(Bound::unbind))?
        .bind(py))
}

fn class<'py>(py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyAny>> {
    errors_module(py)?.getattr(name)
}

/// The Python `Diagnostic` value for one finding.
fn diagnostic_value<'py>(py: Python<'py>, diagnostic: &Diagnostic) -> PyResult<Bound<'py, PyAny>> {
    let span = diagnostic
        .span()
        .map(|span| (span.start.byte_offset, span.end));
    class(py, "Diagnostic")?.call1((
        diagnostic.code().to_string(),
        diagnostic.message(),
        diagnostic.remedy(),
        span,
        diagnostic.severity().label(),
    ))
}

/// Builds `name(message, code=…, remedy=…, span=…, findings=…)`.
fn build(
    py: Python<'_>,
    name: &str,
    message: &str,
    primary: Option<&Diagnostic>,
    findings: &[Diagnostic],
) -> PyResult<PyErr> {
    let keywords = PyDict::new(py);
    if let Some(primary) = primary {
        keywords.set_item("code", primary.code().to_string())?;
        keywords.set_item("remedy", primary.remedy())?;
        keywords.set_item(
            "span",
            primary
                .span()
                .map(|span| (span.start.byte_offset, span.end)),
        )?;
    }
    let values = findings
        .iter()
        .map(|finding| diagnostic_value(py, finding))
        .collect::<PyResult<Vec<_>>>()?;
    keywords.set_item("findings", PyTuple::new(py, values)?)?;
    let instance = class(py, name)?.call((message,), Some(&keywords))?;
    Ok(PyErr::from_value(instance))
}

/// Falls back to a plain `RuntimeError` only when `molframe.errors` itself
/// cannot be used, which would be a broken installation.
fn broken(error: &PyErr, message: &str) -> PyErr {
    pyo3::exceptions::PyRuntimeError::new_err(format!("{message} ({error})"))
}

/// The exception a single diagnostic is raised as.
pub(crate) fn from_diagnostic(diagnostic: &Diagnostic) -> PyErr {
    from_findings(std::slice::from_ref(diagnostic), diagnostic.message())
}

/// The exception for a set of findings, led by the first and carrying all.
///
/// `message` is what `str(error)` reads: a caller that has rendered the
/// findings against their source passes that text, otherwise it is the first
/// finding's own message.
pub(crate) fn from_findings(findings: &[Diagnostic], message: &str) -> PyErr {
    let name = findings
        .first()
        .map_or("MolframeError", |first| class_name(first.code()));
    from_findings_as(name, findings, message)
}

/// As [`from_findings`], raised as the named class instead of the one the
/// first finding's code selects.
pub(crate) fn from_findings_as(name: &str, findings: &[Diagnostic], message: &str) -> PyErr {
    Python::attach(|py| {
        build(py, name, message, findings.first(), findings)
            .unwrap_or_else(|error| broken(&error, message))
    })
}

/// The exception for any kernel error that names a diagnostic.
pub(crate) fn kernel<E>(error: E) -> PyErr
where
    Diagnostic: From<E>,
{
    from_diagnostic(&Diagnostic::from(error))
}

fn plain(name: &str, message: &str) -> PyErr {
    Python::attach(|py| {
        build(py, name, message, None, &[]).unwrap_or_else(|error| broken(&error, message))
    })
}

/// A failure whose kernel names no diagnostic code yet.
///
/// It is still a `MolframeValueError`, so callers catch it with the rest, but
/// its `code` is `None`: classifying it belongs in the kernel's own error type
/// (see `diagnostic_from!`), not here.
pub(crate) fn failure(error: impl std::fmt::Display) -> PyErr {
    plain("MolframeValueError", &error.to_string())
}

/// An argument or input value that is not acceptable.
pub(crate) fn value(message: impl AsRef<str>) -> PyErr {
    plain("MolframeValueError", message.as_ref())
}

/// A name or label that does not exist.
pub(crate) fn key(message: impl AsRef<str>) -> PyErr {
    plain("MolframeKeyError", message.as_ref())
}

/// A position out of range.
pub(crate) fn index(message: impl AsRef<str>) -> PyErr {
    plain("MolframeIndexError", message.as_ref())
}

/// An argument of the wrong type.
pub(crate) fn type_error(message: impl AsRef<str>) -> PyErr {
    plain("MolframeTypeError", message.as_ref())
}

/// A broken internal invariant.
pub(crate) fn internal(message: impl AsRef<str>) -> PyErr {
    plain("InternalError", message.as_ref())
}

/// Raises a non-fatal diagnostic as a warning of `category`.
///
/// The warning is an instance, so `warnings.filterwarnings` can match its
/// class and its message and the handler can read `code` and `remedy`.
pub(crate) fn warn(
    py: Python<'_>,
    category: &str,
    message: &str,
    diagnostic: &Diagnostic,
) -> PyResult<()> {
    let keywords = PyDict::new(py);
    keywords.set_item("code", diagnostic.code().to_string())?;
    keywords.set_item("remedy", diagnostic.remedy())?;
    let warning = class(py, category)?.call((message,), Some(&keywords))?;
    let options = PyDict::new(py);
    options.set_item("stacklevel", 2)?;
    py.import("warnings")?
        .call_method("warn", (warning,), Some(&options))?;
    Ok(())
}

/// Makes `molframe.errors` importable in a test process that embeds Python.
///
/// The extension module is not importable there, so the package is stood in by
/// an empty module whose search path is the source tree; `molframe.errors` is
/// pure Python and needs nothing else.
#[cfg(test)]
pub(crate) fn install_package_for_tests(py: Python<'_>) {
    let source = format!(
        "import sys, types\n\
         package = types.ModuleType('molframe')\n\
         package.__path__ = [{:?}]\n\
         sys.modules.setdefault('molframe', package)\n",
        concat!(env!("CARGO_MANIFEST_DIR"), "/../../python/molframe")
    );
    let code = std::ffi::CString::new(source).expect("no interior nul");
    py.run(&code, None, None)
        .expect("the package stand-in installs");
}

#[cfg(test)]
#[path = "error_tests.rs"]
mod tests;
