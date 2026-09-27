//! Query findings as Python sees them.
//!
//! A query that fails, or that is valid but probably not what was meant, is
//! reported through the one renderer every other surface uses: the query is
//! quoted with the offending token underlined, followed by the details the
//! finding carries (the unknown keyword, the data a keyword needs) and its
//! remedy. A bare "unknown selection keyword" names neither the word nor the
//! fix.

use pyo3::prelude::*;
use std::ffi::CString;

pyo3::create_exception!(
    molframe,
    QueryError,
    pyo3::exceptions::PyValueError,
    "A selection query that cannot be compiled or evaluated."
);

pyo3::create_exception!(
    molframe,
    QueryWarning,
    pyo3::exceptions::PyUserWarning,
    "A selection query that is valid but probably not what was meant."
);

/// Every finding, rendered against the query text it was raised for.
fn render(findings: &[molframe::Diagnostic], source: &str) -> String {
    findings
        .iter()
        .map(|finding| {
            molframe::Rendered::new(finding)
                .with_source(source.as_bytes())
                .to_string()
        })
        .collect::<Vec<_>>()
        .join("\n\n")
}

/// The error a failed query raises.
pub(crate) fn query_error(findings: &molframe::Findings, source: &str) -> PyErr {
    QueryError::new_err(render(findings.as_slice(), source))
}

/// Raises each warning of an evaluated query as a `QueryWarning`, which
/// `warnings.filterwarnings` can silence or promote to an error.
pub(crate) fn warn(
    py: Python<'_>,
    warnings: &[molframe::Diagnostic],
    source: &str,
) -> PyResult<()> {
    let category = py.get_type::<QueryWarning>();
    for warning in warnings {
        let message = CString::new(render(std::slice::from_ref(warning), source))
            .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))?;
        PyErr::warn(py, &category, &message, 2)?;
    }
    Ok(())
}
