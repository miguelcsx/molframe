//! What a structure carries beyond its hierarchy: entry metadata, per-atom
//! annotations and the bulk per-atom columns.

use crate::bindings::PyStructure;
use molframe::{AtomAnnotation, AtomAnnotations};
use numpy::{PyArray1, ToPyArray};
use pyo3::prelude::*;

/// Entry-level facts the file states about the structure.
///
/// A field the file does not state is `None`, never a guess.
#[derive(Clone, Debug)]
#[pyclass(
    name = "EntryMetadata",
    frozen,
    skip_from_py_object,
    module = "molframe"
)]
pub(crate) struct PyEntryMetadata {
    pub(crate) parent: PyStructure,
}

#[pymethods]
impl PyEntryMetadata {
    /// The entry identifier, as deposited.
    #[getter]
    fn id(&self) -> Option<String> {
        self.parent
            .inner
            .metadata()
            .id
            .as_deref()
            .map(str::to_owned)
    }

    #[getter]
    fn title(&self) -> Option<String> {
        self.parent
            .inner
            .metadata()
            .title
            .as_deref()
            .map(str::to_owned)
    }

    /// The experimental method.
    #[getter]
    fn method(&self) -> Option<String> {
        self.parent
            .inner
            .metadata()
            .method
            .as_deref()
            .map(str::to_owned)
    }

    /// Resolution in ångström, where the method reports one.
    #[getter]
    fn resolution(&self) -> Option<f32> {
        self.parent.inner.metadata().resolution
    }

    /// The space-group symbol as the file wrote it.
    #[getter]
    fn space_group(&self) -> Option<String> {
        self.parent
            .inner
            .metadata()
            .space_group
            .as_deref()
            .map(str::to_owned)
    }

    fn __repr__(&self) -> String {
        let text = |value: Option<String>| match value {
            Some(value) => format!("{value:?}"),
            None => "None".to_owned(),
        };
        format!(
            "EntryMetadata(id={}, method={}, resolution={})",
            text(self.id()),
            text(self.method()),
            match self.resolution() {
                Some(value) => value.to_string(),
                None => "None".to_owned(),
            }
        )
    }
}

/// One per-atom column: its values, and which atoms record one.
///
/// `values` is a `NumPy` array (`bool`, `int64` or `float64`) or, for text
/// columns, a list. Where `present` is false the value at that atom is a
/// placeholder and must not be read as data.
#[derive(Debug)]
#[pyclass(name = "Annotation", frozen, skip_from_py_object, module = "molframe")]
pub(crate) struct PyAnnotation {
    kind: &'static str,
    values: Py<PyAny>,
    present: Py<PyArray1<bool>>,
}

#[pymethods]
impl PyAnnotation {
    /// `"boolean"`, `"integer"`, `"real"` or `"symbol"`.
    #[getter]
    const fn kind(&self) -> &'static str {
        self.kind
    }

    #[getter]
    fn values(&self, py: Python<'_>) -> Py<PyAny> {
        self.values.clone_ref(py)
    }

    #[getter]
    fn present(&self, py: Python<'_>) -> Py<PyArray1<bool>> {
        self.present.clone_ref(py)
    }

    fn __repr__(&self) -> String {
        format!("Annotation(kind={})", self.kind)
    }
}

/// Custom per-atom columns by name, in name order.
#[derive(Clone, Debug)]
#[pyclass(name = "Annotations", frozen, skip_from_py_object, module = "molframe")]
pub(crate) struct PyAnnotations {
    pub(crate) parent: PyStructure,
}

impl PyAnnotations {
    fn table(&self) -> &AtomAnnotations {
        self.parent.inner.annotations()
    }
}

fn column(
    py: Python<'_>,
    structure: &PyStructure,
    annotation: &AtomAnnotation,
) -> PyResult<PyAnnotation> {
    macro_rules! typed {
        ($column:expr, $kind:literal) => {{
            let column = $column;
            let present: Vec<bool> = (0..column.len())
                .map(|atom| column.presence(atom).is_present())
                .collect();
            (
                $kind,
                column.values().to_pyarray(py).into_any().unbind(),
                present,
            )
        }};
    }
    let (kind, values, present) = match annotation {
        AtomAnnotation::Boolean(column) => typed!(column, "boolean"),
        AtomAnnotation::Integer(column) => typed!(column, "integer"),
        AtomAnnotation::Real(column) => typed!(column, "real"),
        AtomAnnotation::Symbol(column) => {
            let present: Vec<bool> = (0..column.len())
                .map(|atom| column.presence(atom).is_present())
                .collect();
            let text: Vec<Option<String>> = column
                .values()
                .iter()
                .zip(&present)
                .map(|(symbol, present)| {
                    if *present {
                        structure.inner.engine().resolve(*symbol).map(str::to_owned)
                    } else {
                        None
                    }
                })
                .collect();
            (
                "symbol",
                pyo3::types::PyList::new(py, text)?.into_any().unbind(),
                present,
            )
        }
        _ => {
            return Err(crate::error::type_error(
                "this annotation has a type the Python surface does not carry",
            ));
        }
    };
    Ok(PyAnnotation {
        kind,
        values,
        present: present.to_pyarray(py).unbind(),
    })
}

#[pymethods]
impl PyAnnotations {
    fn __len__(&self) -> usize {
        self.table().len()
    }

    /// The column names, in order.
    #[getter]
    fn names(&self) -> Vec<String> {
        self.table()
            .iter()
            .map(|(name, _)| name.to_owned())
            .collect()
    }

    fn __contains__(&self, name: &str) -> bool {
        self.table().get(name).is_some()
    }

    fn __getitem__(&self, py: Python<'_>, name: &str) -> PyResult<PyAnnotation> {
        let annotation = self
            .table()
            .get(name)
            .ok_or_else(|| crate::error::key(name))?;
        column(py, &self.parent, annotation)
    }

    fn __repr__(&self) -> String {
        format!("Annotations({})", self.names().join(", "))
    }
}

/// The atomic number of every atom.
pub(crate) fn elements(structure: &molframe::Structure) -> Vec<u8> {
    structure
        .engine()
        .data()
        .atoms()
        .map(|atom| match atom.element() {
            Some(element) => element.atomic_number(),
            None => 0,
        })
        .collect()
}

/// The occupancy of every atom, `NaN` where the file records none.
pub(crate) fn occupancies(structure: &molframe::Structure) -> Vec<f32> {
    structure
        .engine()
        .data()
        .atoms()
        .map(|atom| match atom.occupancy() {
            Some(value) => value,
            None => f32::NAN,
        })
        .collect()
}

/// The temperature factor of every atom, `NaN` where the file records none.
pub(crate) fn b_factors(structure: &molframe::Structure) -> Vec<f32> {
    structure
        .engine()
        .data()
        .atoms()
        .map(|atom| match atom.b_factor() {
            Some(value) => value,
            None => f32::NAN,
        })
        .collect()
}

/// The anisotropic displacement table, or `None` when the file never said.
pub(crate) fn anisotropy(
    py: Python<'_>,
    structure: &molframe::Structure,
) -> Option<crate::table::PyTable> {
    let table = structure.anisotropy();
    if !table.is_available() {
        return None;
    }
    let rows: Vec<_> = table.iter().collect();
    let atoms: Vec<u32> = rows.iter().map(|row| row.atom.get()).collect();
    let component =
        |position: usize| -> Vec<f32> { rows.iter().map(|row| row.u[position]).collect() };
    Some(
        crate::table::TableBuilder::new(py, rows.len())
            .indices("atom", &atoms)
            .single("u11", &component(0))
            .single("u22", &component(1))
            .single("u33", &component(2))
            .single("u12", &component(3))
            .single("u13", &component(4))
            .single("u23", &component(5))
            .finish(),
    )
}
