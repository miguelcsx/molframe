//! Conversion of Rust diagnostics into Python exceptions.

use pdbiox::{Class, Code, Diagnostic};
use pyo3::create_exception;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;
use pyo3::sync::PyOnceLock;
use pyo3::types::{PyDict, PyTuple, PyType};

create_exception!(_native, PdbioxError, PyException);
create_exception!(_native, ParseError, PdbioxError);
create_exception!(_native, SchemaError, PdbioxError);
create_exception!(_native, ConsistencyError, PdbioxError);
create_exception!(_native, ConversionError, PdbioxError);
create_exception!(_native, GeometryError, PdbioxError);
create_exception!(_native, PolicyError, PdbioxError);

static INDEX_ERROR: PyOnceLock<Py<PyType>> = PyOnceLock::new();
static KEY_ERROR: PyOnceLock<Py<PyType>> = PyOnceLock::new();

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("PdbioxError", module.py().get_type::<PdbioxError>())?;
    module.add("ParseError", module.py().get_type::<ParseError>())?;
    module.add("SchemaError", module.py().get_type::<SchemaError>())?;
    module.add(
        "ConsistencyError",
        module.py().get_type::<ConsistencyError>(),
    )?;
    module.add("ConversionError", module.py().get_type::<ConversionError>())?;
    module.add("GeometryError", module.py().get_type::<GeometryError>())?;
    module.add("PolicyError", module.py().get_type::<PolicyError>())?;
    let index_error = INDEX_ERROR.get_or_try_init(module.py(), || {
        mixed_exception(module.py(), "PdbioxIndexError", "IndexError")
    })?;
    let key_error = KEY_ERROR.get_or_try_init(module.py(), || {
        mixed_exception(module.py(), "PdbioxKeyError", "KeyError")
    })?;
    module.add("PdbioxIndexError", index_error.bind(module.py()))?;
    module.add("PdbioxKeyError", key_error.bind(module.py()))?;
    Ok(())
}

fn mixed_exception(py: Python<'_>, name: &str, standard: &str) -> PyResult<Py<PyType>> {
    let builtins = py.import("builtins")?;
    let standard = builtins.getattr(standard)?.cast_into::<PyType>()?;
    let bases = PyTuple::new(
        py,
        [py.get_type::<PdbioxError>().into_any(), standard.into_any()],
    )?;
    let namespace = PyDict::new(py);
    Ok(builtins
        .getattr("type")?
        .call1((name, bases, namespace))?
        .cast_into::<PyType>()
        .map(Bound::unbind)?)
}

pub(crate) fn index_error(py: Python<'_>, index: isize) -> PyErr {
    instantiate_standard(py, &INDEX_ERROR, index)
}

pub(crate) fn key_error(py: Python<'_>, key: &str) -> PyErr {
    instantiate_standard(py, &KEY_ERROR, key)
}

fn instantiate_standard(
    py: Python<'_>,
    kind: &PyOnceLock<Py<PyType>>,
    value: impl for<'a> IntoPyObject<'a>,
) -> PyErr {
    let Some(kind) = kind.get(py) else {
        return PdbioxError::new_err(Code::E9001.cause());
    };
    match kind.bind(py).call1((value,)) {
        Ok(instance) => PyErr::from_value(instance),
        Err(error) => error,
    }
}

pub(crate) fn read_error(py: Python<'_>, findings: &[Diagnostic]) -> PyErr {
    let internal_invariant = Diagnostic::new(Code::E9001);
    let finding = match findings.first() {
        Some(finding) => finding,
        None => &internal_invariant,
    };
    let error = match finding.code().class() {
        Class::Syntax => ParseError::new_err(finding.message().to_owned()),
        Class::Schema => SchemaError::new_err(finding.message().to_owned()),
        Class::Consistency => ConsistencyError::new_err(finding.message().to_owned()),
        Class::Conversion => ConversionError::new_err(finding.message().to_owned()),
        Class::Geometry => GeometryError::new_err(finding.message().to_owned()),
        Class::Policy => PolicyError::new_err(finding.message().to_owned()),
        _ => PdbioxError::new_err(finding.message().to_owned()),
    };
    attach_diagnostic(py, error, finding)
}

pub(crate) fn cif_write_error(error: &pdbiox::CifWriteError) -> PyErr {
    ConversionError::new_err(error.to_string())
}

fn attach_diagnostic(py: Python<'_>, error: PyErr, finding: &Diagnostic) -> PyErr {
    let value = error.value(py);
    let attributes = [
        ("code", finding.code().to_string()),
        ("message", finding.message().to_owned()),
        ("remedy", finding.remedy().to_owned()),
    ];
    for (name, attribute) in attributes {
        if let Err(set_error) = value.setattr(name, attribute) {
            return set_error;
        }
    }
    if let Some(span) = finding.span()
        && let Err(set_error) = value.setattr(
            "span",
            (
                span.start.byte_offset,
                span.end,
                span.start.line,
                span.start.column,
            ),
        )
    {
        return set_error;
    }
    error
}
