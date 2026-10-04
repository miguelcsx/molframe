//! Arrow export for named numeric columns of one length.
//!
//! This is the shape of every tabular analysis result: equal-length typed
//! columns with names. Export copies the values once, which the schema says.

use super::extension::{ExportCost, field};
use super::stream::{ArrowStream, ArrowTableExport};
use arrow::array::{
    ArrayRef, Float32Array, Float64Array, Int32Array, Int64Array, UInt8Array, UInt32Array,
};
use arrow::datatypes::{DataType, Schema, SchemaRef};
use arrow::error::{ArrowError, Result};
use arrow::record_batch::RecordBatch;
use std::collections::BTreeSet;
use std::ops::Range;
use std::sync::Arc;

/// One typed column of values.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum Column {
    /// Unsigned bytes, for example secondary-structure codes.
    U8(Vec<u8>),
    /// Unsigned 32-bit integers, for example atom indices.
    U32(Vec<u32>),
    /// Signed 32-bit integers.
    I32(Vec<i32>),
    /// Signed 64-bit integers.
    I64(Vec<i64>),
    /// Single-precision floats.
    F32(Vec<f32>),
    /// Double-precision floats.
    F64(Vec<f64>),
}

impl Column {
    /// Number of values.
    #[must_use]
    pub fn len(&self) -> usize {
        match self {
            Self::U8(values) => values.len(),
            Self::U32(values) => values.len(),
            Self::I32(values) => values.len(),
            Self::I64(values) => values.len(),
            Self::F32(values) => values.len(),
            Self::F64(values) => values.len(),
        }
    }

    /// Whether the column holds no values.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn data_type(&self) -> DataType {
        match self {
            Self::U8(_) => DataType::UInt8,
            Self::U32(_) => DataType::UInt32,
            Self::I32(_) => DataType::Int32,
            Self::I64(_) => DataType::Int64,
            Self::F32(_) => DataType::Float32,
            Self::F64(_) => DataType::Float64,
        }
    }

    fn array(&self) -> ArrayRef {
        match self {
            Self::U8(values) => Arc::new(UInt8Array::from(values.clone())),
            Self::U32(values) => Arc::new(UInt32Array::from(values.clone())),
            Self::I32(values) => Arc::new(Int32Array::from(values.clone())),
            Self::I64(values) => Arc::new(Int64Array::from(values.clone())),
            Self::F32(values) => Arc::new(Float32Array::from(values.clone())),
            Self::F64(values) => Arc::new(Float64Array::from(values.clone())),
        }
    }
}

/// Columns that cannot form a table.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum ColumnTableError {
    /// Two columns share a name.
    #[error("duplicate column name: {name}")]
    DuplicateName {
        /// The repeated name.
        name: String,
    },
    /// A column is not as long as the first.
    #[error("column {name} has {found} values, expected {expected}")]
    LengthMismatch {
        /// The offending column.
        name: String,
        /// Length of the first column.
        expected: usize,
        /// Length of this column.
        found: usize,
    },
}

/// Named columns of one length, exportable as an Arrow stream.
#[derive(Clone, Debug)]
pub struct ColumnTable {
    columns: Arc<[(String, Column)]>,
    length: usize,
    schema: SchemaRef,
}

impl ColumnTable {
    /// Binds columns after checking that names are unique and lengths agree.
    ///
    /// # Errors
    ///
    /// Returns [`ColumnTableError`] for a repeated name or a ragged column.
    pub fn new(columns: Vec<(String, Column)>) -> std::result::Result<Self, ColumnTableError> {
        let length = columns.first().map_or(0, |(_, column)| column.len());
        let mut names = BTreeSet::new();
        for (name, column) in &columns {
            if !names.insert(name.as_str()) {
                return Err(ColumnTableError::DuplicateName { name: name.clone() });
            }
            if column.len() != length {
                return Err(ColumnTableError::LengthMismatch {
                    name: name.clone(),
                    expected: length,
                    found: column.len(),
                });
            }
        }
        let fields = columns
            .iter()
            .map(|(name, column)| field(name, column.data_type(), false, None, ExportCost::Copy))
            .collect::<Vec<_>>();
        Ok(Self {
            columns: columns.into(),
            length,
            schema: Arc::new(Schema::new(fields)),
        })
    }

    /// Number of rows.
    #[must_use]
    pub fn len(&self) -> usize {
        self.length
    }

    /// Whether the table has no rows.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.length == 0
    }

    /// Arrow schema, with the copy cost of every column.
    #[must_use]
    pub fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }

    /// Creates a lazy Arrow C stream over the columns.
    ///
    /// # Errors
    ///
    /// Reserved for failures constructing the stream adapter.
    pub fn arrow_stream(&self) -> Result<ArrowStream> {
        <Self as ArrowTableExport>::arrow_stream(self)
    }

    /// Materialises the columns as record batches.
    ///
    /// # Errors
    ///
    /// Returns an Arrow error if the columns violate the declared schema.
    pub fn record_batches(&self) -> Result<Vec<RecordBatch>> {
        <Self as ArrowTableExport>::record_batches(self)
    }
}

impl ArrowTableExport for ColumnTable {
    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }

    fn batch_count(&self) -> usize {
        usize::from(self.length > 0)
    }

    fn batch(&self, index: usize) -> Result<RecordBatch> {
        if index != 0 || self.length == 0 {
            return Err(ArrowError::InvalidArgumentError(
                "column table batch index is out of bounds".to_owned(),
            ));
        }
        let arrays = self
            .columns
            .iter()
            .map(|(_, column)| column.array())
            .collect();
        RecordBatch::try_new(self.schema.clone(), arrays)
    }

    fn row_count(&self) -> usize {
        self.length
    }

    fn batch_rows(&self, _index: usize) -> Range<usize> {
        0..self.length
    }
}

#[cfg(test)]
#[path = "columns_tests.rs"]
mod tests;
