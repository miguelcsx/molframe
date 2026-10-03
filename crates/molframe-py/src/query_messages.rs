//! Query findings as Python sees them.
//!
//! A query that fails, or that is valid but probably not what was meant, is
//! reported through the one renderer every other surface uses: the query is
//! quoted with the offending token underlined, followed by the details the
//! finding carries (the unknown keyword, the data a keyword needs) and its
//! remedy. A bare "unknown selection keyword" names neither the word nor the
//! fix.

use pyo3::prelude::*;

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

/// The error a failed query raises: a `QueryError` whose text quotes the query
/// and which carries every finding.
pub(crate) fn query_error(findings: &molframe::Findings, source: &str) -> PyErr {
    crate::error::from_findings_as(
        "QueryError",
        findings.as_slice(),
        &render(findings.as_slice(), source),
    )
}

/// Raises each warning of an evaluated query as a `QueryWarning`, which
/// `warnings.filterwarnings` can silence or promote to an error.
pub(crate) fn warn(
    py: Python<'_>,
    warnings: &[molframe::Diagnostic],
    source: &str,
) -> PyResult<()> {
    for warning in warnings {
        let message = render(std::slice::from_ref(warning), source);
        crate::error::warn(py, "QueryWarning", &message, warning)?;
    }
    Ok(())
}
