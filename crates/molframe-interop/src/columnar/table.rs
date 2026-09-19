//! Shared plumbing for the derived residue, chain and bond tables.
//!
//! Each derived table batches at a fixed [`TABLE_BATCH_ROWS`] row count and
//! materialises one batch per pull, so cost is `O(rows)` per batch and
//! `O(1)` to construct.

use super::extension::{ExportCost, field};

pub(crate) const TABLE_BATCH_ROWS: usize = 65_536;

macro_rules! table_type {
    ($name:ident, $schema:ident, $count:ident, $batch:ident, $description:literal) => {
        #[doc = $description]
        #[doc = "\n\nCloning is `O(1)`: the immutable snapshot is shared through [`Arc`]."]
        #[derive(Clone, Debug)]
        pub struct $name {
            structure: ::std::sync::Arc<::molframe_core::Structure>,
            schema: ::std::sync::Arc<::arrow::datatypes::Schema>,
        }

        impl $name {
            /// Binds an immutable structure snapshot to this Arrow table.
            #[must_use]
            pub fn new(structure: &::molframe_core::Structure) -> Self {
                Self {
                    structure: ::std::sync::Arc::new(structure.clone()),
                    schema: ::std::sync::Arc::new($schema()),
                }
            }

            /// Arrow schema with molframe extension and copy-cost metadata.
            #[must_use]
            pub fn schema(&self) -> ::arrow::datatypes::SchemaRef {
                self.schema.clone()
            }

            /// Materialises the derived table as record batches.
            ///
            /// # Errors
            ///
            /// Returns an Arrow error if columns violate the declared schema.
            pub fn record_batches(
                &self,
            ) -> ::arrow::error::Result<Vec<::arrow::record_batch::RecordBatch>> {
                <$name as crate::columnar::stream::ArrowTableExport>::record_batches(self)
            }

            /// Creates a consumable Arrow C Stream Interface value.
            ///
            /// The derived table is converted only when the consumer pulls its
            /// batch. A conversion error is reported by that pull.
            ///
            /// # Errors
            ///
            /// Reserved for failures constructing the stream adapter.
            pub fn arrow_stream(
                &self,
            ) -> ::arrow::error::Result<crate::columnar::stream::ArrowStream> {
                <$name as crate::columnar::stream::ArrowTableExport>::arrow_stream(self)
            }
        }

        impl crate::columnar::stream::ArrowTableExport for $name {
            fn schema(&self) -> ::arrow::datatypes::SchemaRef {
                self.schema.clone()
            }

            fn batch_count(&self) -> usize {
                $count(&self.structure).div_ceil(crate::columnar::table::TABLE_BATCH_ROWS)
            }

            fn batch(
                &self,
                index: usize,
            ) -> ::arrow::error::Result<::arrow::record_batch::RecordBatch> {
                let count = $count(&self.structure);
                let start = index
                    .checked_mul(crate::columnar::table::TABLE_BATCH_ROWS)
                    .ok_or_else(|| {
                        ::arrow::error::ArrowError::InvalidArgumentError(
                            "table batch index overflows".to_owned(),
                        )
                    })?;
                if start >= count {
                    return Err(::arrow::error::ArrowError::InvalidArgumentError(
                        "table batch index is out of bounds".to_owned(),
                    ));
                }
                $batch(
                    &self.structure,
                    self.schema.clone(),
                    start
                        ..start
                            .saturating_add(crate::columnar::table::TABLE_BATCH_ROWS)
                            .min(count),
                )
            }
        }
    };
}

pub(crate) use table_type;

/// Returns an Arrow field marked as decoded in one pass.
pub(crate) fn decoded(
    name: &str,
    data_type: ::arrow::datatypes::DataType,
    extension: &str,
    nullable: bool,
) -> ::arrow::datatypes::Field {
    field(
        name,
        data_type,
        nullable,
        Some(extension),
        ExportCost::Decode,
    )
}

/// Width of a topology range, checked against a reversed range.
pub(crate) fn range_width(range: std::ops::Range<u32>) -> ::arrow::error::Result<u32> {
    range.end.checked_sub(range.start).ok_or_else(|| {
        ::arrow::error::ArrowError::InvalidArgumentError(
            "topology range ends before it starts".to_owned(),
        )
    })
}
