//! Read-only Python projections for lossless mmCIF documents.

use crate::cif_rows::PyCifRows;
use crate::errors::read_error;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::path::PathBuf;

#[pyclass(name = "Quoting", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyCifQuoting {
    Bare,
    Single,
    Double,
    Text,
}

impl From<PyCifQuoting> for pdbiox::cif::lexer::Quoting {
    fn from(value: PyCifQuoting) -> Self {
        match value {
            PyCifQuoting::Bare => Self::Bare,
            PyCifQuoting::Single => Self::Single,
            PyCifQuoting::Double => Self::Double,
            PyCifQuoting::Text => Self::Text,
        }
    }
}

impl From<pdbiox::cif::lexer::Quoting> for PyCifQuoting {
    fn from(value: pdbiox::cif::lexer::Quoting) -> Self {
        match value {
            pdbiox::cif::lexer::Quoting::Bare => Self::Bare,
            pdbiox::cif::lexer::Quoting::Single => Self::Single,
            pdbiox::cif::lexer::Quoting::Double => Self::Double,
            pdbiox::cif::lexer::Quoting::Text => Self::Text,
        }
    }
}

#[pyclass(name = "CifValue", from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyCifValue {
    pub(crate) inner: pdbiox::CifValue,
}

#[pymethods]
impl PyCifValue {
    #[new]
    fn new(text: &str, quoting: PyCifQuoting) -> Self {
        Self {
            inner: pdbiox::CifValue::parse(text, quoting.into()),
        }
    }

    #[staticmethod]
    fn parse(text: &str, quoting: PyCifQuoting) -> Self {
        Self::new(text, quoting)
    }

    #[staticmethod]
    fn inapplicable() -> Self {
        Self {
            inner: pdbiox::CifValue::Inapplicable,
        }
    }

    #[staticmethod]
    fn unknown() -> Self {
        Self {
            inner: pdbiox::CifValue::Unknown,
        }
    }

    #[staticmethod]
    fn text_value(value: String) -> Self {
        Self {
            inner: pdbiox::CifValue::Text(value.into()),
        }
    }

    #[staticmethod]
    fn integer_value(value: i64) -> Self {
        Self {
            inner: pdbiox::CifValue::Integer(value),
        }
    }

    #[staticmethod]
    fn float_value(value: f64) -> PyResult<Self> {
        if !value.is_finite() {
            return Err(PyValueError::new_err("CIF float values must be finite"));
        }
        Ok(Self {
            inner: pdbiox::CifValue::Float(value),
        })
    }

    #[getter]
    fn kind(&self) -> &'static str {
        match &self.inner {
            pdbiox::CifValue::Inapplicable => "inapplicable",
            pdbiox::CifValue::Unknown => "unknown",
            pdbiox::CifValue::Text(_) => "text",
            pdbiox::CifValue::Integer(_) => "integer",
            pdbiox::CifValue::Float(_) => "float",
        }
    }

    #[getter]
    fn text(&self) -> Option<String> {
        self.inner.as_str().map(str::to_owned)
    }

    fn as_str(&self) -> Option<String> {
        self.text()
    }

    #[getter]
    fn integer(&self) -> Option<i64> {
        self.inner.as_integer()
    }

    #[getter]
    fn number(&self) -> Option<f64> {
        self.inner.as_float()
    }

    fn as_integer(&self) -> Option<i64> {
        self.inner.as_integer()
    }

    fn as_float(&self) -> Option<f64> {
        self.inner.as_float()
    }

    fn as_identifier(&self) -> Option<String> {
        self.inner.as_identifier().map(std::borrow::Cow::into_owned)
    }

    fn is_recorded(&self) -> bool {
        self.inner.is_recorded()
    }
}

#[pyclass(name = "Column", from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyCifColumn {
    inner: pdbiox::Column,
}

#[pymethods]
impl PyCifColumn {
    #[new]
    #[pyo3(signature = (values=None))]
    fn new(values: Option<Vec<(PyCifValue, PyCifQuoting)>>) -> Self {
        let mut column = Self {
            inner: pdbiox::Column::default(),
        };
        if let Some(values) = values {
            for (value, quoting) in values {
                column.inner.push(value.inner, quoting.into());
            }
        }
        column
    }

    fn push(&mut self, value: PyCifValue, quoting: PyCifQuoting) {
        self.inner.push(value.inner, quoting.into());
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __getitem__(&self, index: isize) -> Option<PyCifValue> {
        let index = usize::try_from(index).ok()?;
        self.inner.get(index).cloned().map(Into::into)
    }

    fn get(&self, row: usize) -> Option<PyCifValue> {
        self.inner.get(row).cloned().map(Into::into)
    }

    fn quoting(&self, row: usize) -> Option<PyCifQuoting> {
        self.inner.quoting(row).map(Into::into)
    }

    #[getter]
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    #[getter]
    fn values(&self) -> Vec<PyCifValue> {
        self.inner.iter().cloned().map(Into::into).collect()
    }
}

#[pyclass(name = "Category", from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyCifCategory {
    pub(crate) inner: pdbiox::Category,
}

#[pymethods]
impl PyCifCategory {
    #[new]
    #[pyo3(signature = (name, span=None))]
    fn new(name: String, span: Option<(u32, u32, u32, u32)>) -> Self {
        Self {
            inner: pdbiox::Category::new(name, span.map_or_else(default_span, span_value)),
        }
    }

    #[getter]
    fn name(&self) -> &str {
        self.inner.name()
    }

    #[getter]
    fn row_count(&self) -> usize {
        self.inner.row_count()
    }

    #[getter]
    fn len(&self) -> usize {
        self.inner.len()
    }

    #[getter]
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    #[getter]
    fn span(&self) -> (u32, u32, u32, u32) {
        span_tuple(self.inner.span())
    }

    #[getter]
    fn items(&self) -> Vec<String> {
        self.inner.items().map(str::to_owned).collect()
    }

    fn column(&self, item: &str) -> Option<PyCifColumn> {
        self.inner.column(item).cloned().map(Into::into)
    }

    fn column_mut(&self, item: &str) -> Option<PyCifColumn> {
        self.column(item)
    }

    fn rows(&self) -> PyCifRows {
        PyCifRows::from_category(self.inner.clone())
    }

    fn push(&mut self, item: &str, value: PyCifValue, quoting: PyCifQuoting) {
        self.inner
            .column_mut(item)
            .push(value.inner, quoting.into());
    }

    fn value(&self, item: &str, row: usize) -> Option<PyCifValue> {
        self.inner.value(item, row).cloned().map(Into::into)
    }

    fn text(&self, item: &str, row: usize) -> Option<String> {
        self.inner.text(item, row).map(str::to_owned)
    }

    fn identifier(&self, item: &str, row: usize) -> Option<String> {
        self.inner
            .identifier(item, row)
            .map(std::borrow::Cow::into_owned)
    }
}

#[pyclass(name = "DataBlock", from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyCifDataBlock {
    inner: pdbiox::DataBlock,
}

#[pymethods]
impl PyCifDataBlock {
    #[new]
    fn new(name: String) -> Self {
        Self {
            inner: pdbiox::DataBlock::new(name),
        }
    }

    #[getter]
    fn name(&self) -> &str {
        self.inner.name()
    }

    #[getter]
    fn categories(&self) -> Vec<PyCifCategory> {
        self.inner.categories().cloned().map(Into::into).collect()
    }

    #[getter]
    fn len(&self) -> usize {
        self.inner.len()
    }

    #[getter]
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    fn category(&self, name: &str) -> Option<PyCifCategory> {
        self.inner.category(name).cloned().map(Into::into)
    }

    fn category_mut(&self, name: &str) -> Option<PyCifCategory> {
        self.category(name)
    }

    fn add_category(&mut self, category: PyCifCategory) {
        let name = category.inner.name().to_owned();
        let target = self.inner.category_mut(&name, category.inner.span());
        for item in category.inner.items() {
            if let Some(column) = category.inner.column(item) {
                for row in 0..column.len() {
                    if let Some(value) = column.get(row).cloned() {
                        let quoting = match column.quoting(row) {
                            Some(quoting) => quoting,
                            None => pdbiox::cif::lexer::Quoting::Bare,
                        };
                        target.column_mut(item).push(value, quoting);
                    }
                }
            }
        }
    }

    #[pyo3(signature = (category, item, value, quoting, span=None))]
    fn push(
        &mut self,
        category: &str,
        item: &str,
        value: PyCifValue,
        quoting: PyCifQuoting,
        span: Option<(u32, u32, u32, u32)>,
    ) {
        self.inner
            .category_mut(category, span.map_or_else(default_span, span_value))
            .column_mut(item)
            .push(value.inner, quoting.into());
    }
}

#[pyclass(name = "Document", from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyCifDocument {
    pub(crate) inner: pdbiox::Document,
}

#[pymethods]
impl PyCifDocument {
    #[new]
    fn new() -> Self {
        Self {
            inner: pdbiox::Document::new(),
        }
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    #[getter]
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    fn push(&mut self, block: PyCifDataBlock) {
        self.inner.push(block.inner);
    }

    #[getter]
    fn blocks(&self) -> Vec<PyCifDataBlock> {
        self.inner.blocks().cloned().map(Into::into).collect()
    }

    fn first_block(&self) -> Option<PyCifDataBlock> {
        self.inner.first_block().cloned().map(Into::into)
    }

    fn last_block(&self) -> Option<PyCifDataBlock> {
        self.inner.blocks().last().cloned().map(Into::into)
    }

    fn last_block_mut(&self) -> Option<PyCifDataBlock> {
        self.last_block()
    }
}

#[pyfunction]
pub(crate) fn read_document(py: Python<'_>, path: PathBuf) -> PyResult<PyCifDocument> {
    py.detach(move || pdbiox::read_document(path))
        .map(|inner| PyCifDocument { inner })
        .map_err(|findings| read_error(py, &findings))
}

#[pyfunction]
pub(crate) fn write_preserving(document: &PyCifDocument) -> String {
    pdbiox::write_preserving(&document.inner)
}

impl From<pdbiox::CifValue> for PyCifValue {
    fn from(inner: pdbiox::CifValue) -> Self {
        Self { inner }
    }
}

impl From<pdbiox::Document> for PyCifDocument {
    fn from(inner: pdbiox::Document) -> Self {
        Self { inner }
    }
}

impl From<pdbiox::Column> for PyCifColumn {
    fn from(inner: pdbiox::Column) -> Self {
        Self { inner }
    }
}

impl From<pdbiox::Category> for PyCifCategory {
    fn from(inner: pdbiox::Category) -> Self {
        Self { inner }
    }
}

impl From<pdbiox::DataBlock> for PyCifDataBlock {
    fn from(inner: pdbiox::DataBlock) -> Self {
        Self { inner }
    }
}

fn default_span() -> pdbiox::ByteSpan {
    pdbiox::ByteSpan::default()
}

fn span_value(value: (u32, u32, u32, u32)) -> pdbiox::ByteSpan {
    pdbiox::ByteSpan::new(pdbiox::Position::new(value.0, value.1, value.2), value.3)
}

fn span_tuple(value: pdbiox::ByteSpan) -> (u32, u32, u32, u32) {
    (
        value.start.byte_offset,
        value.start.line,
        value.start.column,
        value.end,
    )
}

#[cfg(test)]
#[path = "cif_document_tests.rs"]
mod tests;
