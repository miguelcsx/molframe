//! Owned Python cursor adapting the borrowed CIF `Rows` iterator.

use crate::cif_document::{PyCifCategory, PyCifValue};
use pyo3::prelude::*;

#[pyclass(name = "Rows", from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyCifRows {
    category: molframe::Category,
    row: usize,
}

#[pymethods]
impl PyCifRows {
    #[new]
    fn new(category: PyCifCategory) -> Self {
        Self {
            category: category.inner,
            row: 0,
        }
    }

    #[getter]
    fn category(&self) -> PyCifCategory {
        self.category.clone().into()
    }

    #[getter]
    const fn row(&self) -> usize {
        self.row
    }

    fn advance(&mut self) -> bool {
        self.row = self.row.saturating_add(1);
        self.row < self.category.row_count()
    }

    fn value(&self, item: &str) -> Option<PyCifValue> {
        self.category.value(item, self.row).cloned().map(Into::into)
    }

    fn text(&self, item: &str) -> Option<String> {
        self.category.text(item, self.row).map(str::to_owned)
    }

    fn identifier(&self, item: &str) -> Option<String> {
        self.category
            .identifier(item, self.row)
            .map(std::borrow::Cow::into_owned)
    }

    fn integer(&self, item: &str) -> Option<i64> {
        self.category
            .value(item, self.row)
            .and_then(molframe::CifValue::as_integer)
    }

    fn float(&self, item: &str) -> Option<f64> {
        self.category
            .value(item, self.row)
            .and_then(molframe::CifValue::as_float)
    }

    fn is_recorded(&self, item: &str) -> bool {
        self.category
            .value(item, self.row)
            .is_some_and(molframe::CifValue::is_recorded)
    }
}

impl PyCifRows {
    pub(crate) fn from_category(category: molframe::Category) -> Self {
        Self { category, row: 0 }
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyCifRows>()
}
