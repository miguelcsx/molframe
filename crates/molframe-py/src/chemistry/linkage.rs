//! Mechanical conversion of explicit polymer connectivity policies.

use molframe::chemistry as chem;
use pyo3::prelude::*;

#[derive(Clone)]
#[pyclass(
    name = "PolymerLinkRule",
    frozen,
    skip_from_py_object,
    module = "molframe.chemistry"
)]
struct PyPolymerLinkRule {
    inner: chem::PolymerLinkRule,
}

#[pymethods]
impl PyPolymerLinkRule {
    #[new]
    fn new(
        left_kind: i64,
        left_atom: String,
        right_kind: i64,
        right_atom: String,
    ) -> PyResult<Self> {
        let kind = |code| {
            chem::ComponentKind::from_code(code)
                .ok_or_else(|| crate::error::value("invalid component kind code"))
        };
        Ok(Self {
            inner: chem::PolymerLinkRule::new(
                kind(left_kind)?,
                left_atom,
                kind(right_kind)?,
                right_atom,
            ),
        })
    }

    #[getter]
    fn left_kind(&self) -> i64 {
        self.inner.left_kind.code()
    }
    #[getter]
    fn left_atom(&self) -> &str {
        &self.inner.left_atom
    }
    #[getter]
    fn right_kind(&self) -> i64 {
        self.inner.right_kind.code()
    }
    #[getter]
    fn right_atom(&self) -> &str {
        &self.inner.right_atom
    }
}

/// Explicit distance and attachment rules; annotation validates their chemistry.
#[derive(Clone)]
#[pyclass(
    name = "PolymerLinkPolicy",
    frozen,
    skip_from_py_object,
    module = "molframe.chemistry"
)]
pub(super) struct PyPolymerLinkPolicy {
    pub(super) inner: chem::PolymerLinkPolicy,
}

#[pymethods]
impl PyPolymerLinkPolicy {
    #[new]
    fn new(angstrom: f32, rules: Vec<PyRef<'_, PyPolymerLinkRule>>) -> Self {
        Self {
            inner: chem::PolymerLinkPolicy::explicit(
                angstrom,
                rules
                    .into_iter()
                    .map(|rule| rule.inner.clone())
                    .collect::<Vec<_>>(),
            ),
        }
    }

    #[getter]
    fn angstrom(&self) -> Option<f32> {
        match &self.inner {
            chem::PolymerLinkPolicy::Disabled => None,
            chem::PolymerLinkPolicy::Explicit { angstrom, .. } => Some(*angstrom),
        }
    }

    #[getter]
    fn rules(&self) -> Vec<PyPolymerLinkRule> {
        match &self.inner {
            chem::PolymerLinkPolicy::Disabled => Vec::new(),
            chem::PolymerLinkPolicy::Explicit { rules, .. } => rules
                .iter()
                .map(|rule| PyPolymerLinkRule {
                    inner: rule.clone(),
                })
                .collect(),
        }
    }
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyPolymerLinkRule>()?;
    module.add_class::<PyPolymerLinkPolicy>()
}
