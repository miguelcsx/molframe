//! Complete Python projection of the core diagnostic vocabulary.

use pyo3::prelude::*;

#[pyclass(name = "Kind", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyKind {
    Error,
    Warning,
}

impl From<molframe::Kind> for PyKind {
    fn from(value: molframe::Kind) -> Self {
        match value {
            molframe::Kind::Error => Self::Error,
            molframe::Kind::Warning => Self::Warning,
        }
    }
}

impl From<PyKind> for molframe::Kind {
    fn from(value: PyKind) -> Self {
        match value {
            PyKind::Error => Self::Error,
            PyKind::Warning => Self::Warning,
        }
    }
}

#[pyclass(name = "Severity", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum PySeverity {
    Info,
    Loose,
    Strict,
    Invalidating,
    Breaking,
}

impl From<molframe::Severity> for PySeverity {
    fn from(value: molframe::Severity) -> Self {
        match value {
            molframe::Severity::Info => Self::Info,
            molframe::Severity::Loose => Self::Loose,
            molframe::Severity::Strict => Self::Strict,
            molframe::Severity::Invalidating => Self::Invalidating,
            molframe::Severity::Breaking => Self::Breaking,
        }
    }
}

impl From<PySeverity> for molframe::Severity {
    fn from(value: PySeverity) -> Self {
        match value {
            PySeverity::Info => Self::Info,
            PySeverity::Loose => Self::Loose,
            PySeverity::Strict => Self::Strict,
            PySeverity::Invalidating => Self::Invalidating,
            PySeverity::Breaking => Self::Breaking,
        }
    }
}

#[pymethods]
impl PySeverity {
    fn is_error(&self, strictness: PyStrictness) -> bool {
        molframe::Severity::from(*self).is_error(strictness.into())
    }

    #[getter]
    fn label(&self) -> &'static str {
        molframe::Severity::from(*self).label()
    }
}

#[pyclass(name = "Strictness", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyStrictness {
    Strict,
    Medium,
    Loose,
}

impl From<PyStrictness> for molframe::Strictness {
    fn from(value: PyStrictness) -> Self {
        match value {
            PyStrictness::Strict => Self::Strict,
            PyStrictness::Medium => Self::Medium,
            PyStrictness::Loose => Self::Loose,
        }
    }
}

#[pyclass(name = "Class", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyClass {
    Syntax,
    Schema,
    Consistency,
    Conversion,
    Geometry,
    Policy,
    Resource,
    Internal,
}

impl From<molframe::Class> for PyClass {
    fn from(value: molframe::Class) -> Self {
        match value {
            molframe::Class::Syntax => Self::Syntax,
            molframe::Class::Schema => Self::Schema,
            molframe::Class::Consistency => Self::Consistency,
            molframe::Class::Conversion => Self::Conversion,
            molframe::Class::Geometry => Self::Geometry,
            molframe::Class::Policy => Self::Policy,
            molframe::Class::Resource => Self::Resource,
            _ => Self::Internal,
        }
    }
}

#[pyclass(name = "Code", frozen, eq, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct PyCode(pub(crate) molframe::Code);

#[pymethods]
impl PyCode {
    #[new]
    fn new(kind: PyKind, number: u16) -> Self {
        Self(molframe::Code::new(kind.into(), number))
    }

    #[staticmethod]
    fn from_text(text: &str) -> Option<Self> {
        let text = match text.strip_prefix("MOLFRAME-") {
            Some(value) => value,
            None => text,
        };
        let (kind, number) = text.split_at(1);
        let kind = match kind {
            "E" => PyKind::Error,
            "W" => PyKind::Warning,
            _ => return None,
        };
        let number = number.parse::<u16>().ok()?;
        Some(Self::new(kind, number))
    }

    #[staticmethod]
    fn registered() -> Vec<Self> {
        molframe::Code::registered().map(Self).collect()
    }

    #[getter]
    fn kind(&self) -> PyKind {
        self.0.kind().into()
    }

    #[getter]
    fn number(&self) -> u16 {
        self.0.number()
    }

    #[getter]
    #[pyo3(name = "class")]
    fn class_(&self) -> PyClass {
        self.0.class().into()
    }

    #[getter]
    #[pyo3(name = "class_")]
    fn class_alias(&self) -> PyClass {
        self.0.class().into()
    }

    #[getter]
    fn severity(&self) -> PySeverity {
        self.0.severity().into()
    }

    #[getter]
    fn cause(&self) -> &'static str {
        self.0.cause()
    }

    #[getter]
    fn remedy(&self) -> &'static str {
        self.0.remedy()
    }

    #[getter]
    fn is_registered(&self) -> bool {
        self.0.is_registered()
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("Code('{}')", self.0)
    }
}

#[pyclass(name = "ContextItem", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyContextItem {
    #[pyo3(get)]
    pub(crate) label: String,
    #[pyo3(get)]
    pub(crate) value: String,
}

impl From<molframe::ContextItem> for PyContextItem {
    fn from(value: molframe::ContextItem) -> Self {
        Self {
            label: value.label().to_owned(),
            value: value.value().to_owned(),
        }
    }
}

#[pyclass(name = "Rendered", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyRendered {
    #[pyo3(get)]
    pub(crate) text: String,
}

#[pymethods]
impl PyRendered {
    fn __str__(&self) -> &str {
        &self.text
    }

    fn __repr__(&self) -> String {
        format!("Rendered({:?})", self.text)
    }
}

#[pyclass(name = "Diagnostics", from_py_object)]
#[derive(Clone, Debug, Default)]
pub(crate) struct PyDiagnostics(pub(crate) molframe::Diagnostics);

#[pymethods]
impl PyDiagnostics {
    #[new]
    #[pyo3(signature = (capacity=0))]
    fn new(capacity: usize) -> Self {
        if capacity == 0 {
            Self(molframe::Diagnostics::new())
        } else {
            Self(molframe::Diagnostics::with_capacity(capacity))
        }
    }

    fn push(&mut self, finding: &crate::contract::PyDiagnostic) {
        self.0.push(finding.inner.clone());
    }

    #[getter]
    fn worst(&self) -> Option<PySeverity> {
        self.0.worst().map(Into::into)
    }

    fn has_error(&self, strictness: PyStrictness) -> bool {
        self.0.has_error(strictness.into())
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    fn as_list(&self) -> Vec<crate::contract::PyDiagnostic> {
        self.0.as_slice().iter().cloned().map(Into::into).collect()
    }

    fn finish(&mut self) -> Vec<crate::contract::PyDiagnostic> {
        std::mem::take(&mut self.0)
            .finish()
            .into_iter()
            .map(Into::into)
            .collect()
    }
}

impl PyRendered {
    pub(crate) fn from_diagnostic(
        diagnostic: &molframe::Diagnostic,
        source: Option<&[u8]>,
        origin: Option<&str>,
        color: bool,
    ) -> Self {
        let rendered = match (source, origin) {
            (Some(source), Some(origin)) => molframe::Rendered::new(diagnostic)
                .with_source(source)
                .with_origin(origin)
                .with_color(color)
                .to_string(),
            (Some(source), None) => molframe::Rendered::new(diagnostic)
                .with_source(source)
                .with_color(color)
                .to_string(),
            (None, Some(origin)) => molframe::Rendered::new(diagnostic)
                .with_origin(origin)
                .with_color(color)
                .to_string(),
            (None, None) => molframe::Rendered::new(diagnostic)
                .with_color(color)
                .to_string(),
        };
        Self { text: rendered }
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyKind>()?;
    module.add_class::<PySeverity>()?;
    module.add_class::<PyStrictness>()?;
    module.add_class::<PyClass>()?;
    module.add_class::<PyCode>()?;
    module.add_class::<PyContextItem>()?;
    module.add_class::<PyRendered>()?;
    module.add_class::<PyDiagnostics>()?;
    Ok(())
}
