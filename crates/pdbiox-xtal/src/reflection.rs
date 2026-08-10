//! Shared reflection-table model for structure-factor mmCIF and MTZ.

use pdbiox_core::structure::UnitCell;

use crate::numeric::{f64_to_i32, i64_to_f64};

/// One missing-aware reflection-table value.
#[derive(Clone, Debug, PartialEq)]
pub enum ReflectionValue {
    /// Value was not recorded (`?` in mmCIF or a missing MTZ number).
    Missing,
    /// Value does not apply (`.` in mmCIF).
    Inapplicable,
    /// Exact integer value.
    Integer(i64),
    /// Floating-point measurement.
    Real(f64),
    /// Text or status value.
    Text(Box<str>),
}

impl ReflectionValue {
    /// Returns a numeric value, preserving integers exactly until conversion.
    #[must_use]
    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Self::Integer(value) => Some(i64_to_f64(*value)),
            Self::Real(value) => Some(*value),
            Self::Missing | Self::Inapplicable | Self::Text(_) => None,
        }
    }

    /// Returns an integer only when the stored numeric value is integral and in range.
    #[must_use]
    pub fn as_i32(&self) -> Option<i32> {
        match self {
            Self::Integer(value) => i32::try_from(*value).ok(),
            Self::Real(value)
                if value.is_finite()
                    && value.fract().abs() <= f64::EPSILON
                    && *value >= f64::from(i32::MIN)
                    && *value <= f64::from(i32::MAX) =>
            {
                Some(f64_to_i32(*value))
            }
            Self::Missing | Self::Inapplicable | Self::Real(_) | Self::Text(_) => None,
        }
    }
}

/// Semantic type of one reflection column.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReflectionColumnType {
    /// Miller index.
    MillerIndex,
    /// Structure-factor amplitude.
    Amplitude,
    /// Intensity or squared amplitude.
    Intensity,
    /// Standard uncertainty.
    StandardDeviation,
    /// Phase in degrees.
    Phase,
    /// Free-set or other integer flag.
    Flag,
    /// Generic real-valued column.
    Real,
    /// Text/status column, available in mmCIF but not an MTZ data column.
    Text,
}

/// One named reflection column.
#[derive(Clone, Debug, PartialEq)]
pub struct ReflectionColumn {
    /// Format-native label (`index_h`, `F_meas_au`, `H`, `FP`, ...).
    pub label: Box<str>,
    /// Semantic column type.
    pub column_type: ReflectionColumnType,
    /// One value per reflection row.
    pub values: Vec<ReflectionValue>,
    /// Dataset identifier for multi-dataset MTZ files.
    pub dataset_id: i32,
    /// Original one-character MTZ type code, when the source was MTZ.
    pub mtz_type: Option<char>,
}

/// Metadata for one MTZ/mmCIF diffraction dataset.
#[derive(Clone, Debug, PartialEq)]
pub struct ReflectionDataset {
    /// Stable dataset identifier.
    pub id: i32,
    /// Project name.
    pub project: Box<str>,
    /// Crystal name.
    pub crystal: Box<str>,
    /// Dataset name.
    pub name: Box<str>,
    /// Experimental wavelength in ångströms, when known.
    pub wavelength: Option<f64>,
    /// Dataset-specific cell, when present.
    pub cell: Option<UnitCell>,
}

/// A missing-aware table of reciprocal-space observations.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReflectionTable {
    /// Human-readable title or source block name.
    pub title: Box<str>,
    /// Crystal unit cell.
    pub cell: Option<UnitCell>,
    /// International Tables space-group number.
    pub space_group_number: Option<i32>,
    /// Hermann-Mauguin space-group name.
    pub space_group_name: Option<Box<str>>,
    /// Reflection columns in source order.
    pub columns: Vec<ReflectionColumn>,
    /// Dataset definitions.
    pub datasets: Vec<ReflectionDataset>,
    /// Ordered provenance/history lines.
    pub history: Vec<Box<str>>,
    /// Symmetry operations in International Tables triplet notation.
    pub symmetry_operations: Vec<Box<str>>,
    /// MTZ sort order, using one-based column numbers and zero for unused slots.
    pub sort_order: [i32; 5],
    /// Minimum and maximum reciprocal resolution stored as `1/d²`.
    pub resolution_range: Option<[f64; 2]>,
    /// Explicit MTZ missing-number sentinel; `None` selects IEEE NaN.
    pub missing_value: Option<f32>,
    /// Unrecognised 80-character MTZ main-header records retained in order.
    pub extra_header_records: Vec<Box<str>>,
}

/// Invalid reflection table or required-column request.
#[derive(Clone, Debug, thiserror::Error, PartialEq, Eq)]
#[non_exhaustive]
pub enum ReflectionError {
    /// No reflection category or columns were present.
    #[error("reflection input has no reflection table")]
    MissingTable,
    /// Columns do not all have the same non-zero row count.
    #[error("reflection columns have inconsistent lengths")]
    ColumnLength,
    /// A required Miller-index column is absent or contains a non-integer value.
    #[error("reflection table has invalid or missing Miller indices")]
    MillerIndices,
    /// Crystal metadata is malformed.
    #[error("reflection table has invalid crystal metadata")]
    Metadata,
    /// A numeric value cannot be represented by the target format.
    #[error("reflection table contains an invalid numeric value")]
    InvalidNumber,
    /// An input format feature is unsupported without losing information.
    #[error("unsupported reflection format feature: {0}")]
    Unsupported(Box<str>),
    /// MTZ framing or header records are malformed or inconsistent.
    #[error("invalid or inconsistent MTZ reflection file")]
    InvalidMtz,
    /// MTZ data ends before the declared reflections or headers.
    #[error("truncated MTZ reflection file")]
    TruncatedMtz,
}

impl ReflectionTable {
    /// Number of reflection rows.
    #[must_use]
    pub fn row_count(&self) -> usize {
        self.columns.first().map_or(0, |column| column.values.len())
    }

    /// Validates rectangularity and finite recorded real values.
    ///
    /// # Errors
    ///
    /// Returns an error for an empty/ragged table or non-finite number.
    pub fn validate(&self) -> Result<(), ReflectionError> {
        let rows = self.row_count();
        if rows == 0
            || self
                .columns
                .iter()
                .any(|column| column.values.len() != rows)
        {
            return Err(ReflectionError::ColumnLength);
        }
        if self
            .columns
            .iter()
            .flat_map(|column| &column.values)
            .any(|value| matches!(value, ReflectionValue::Real(number) if !number.is_finite()))
        {
            return Err(ReflectionError::InvalidNumber);
        }
        Ok(())
    }

    /// Finds a column case-insensitively by any of the supplied labels.
    #[must_use]
    pub fn column_any(&self, labels: &[&str]) -> Option<&ReflectionColumn> {
        self.columns.iter().find(|column| {
            labels
                .iter()
                .any(|label| column.label.eq_ignore_ascii_case(label))
        })
    }

    /// Returns Miller indices from mmCIF (`index_h/k/l`) or MTZ (`H/K/L`).
    ///
    /// # Errors
    ///
    /// Returns an error when a required column is absent or non-integral.
    pub fn miller_indices(&self) -> Result<Vec<[i32; 3]>, ReflectionError> {
        let h = self
            .column_any(&["index_h", "H"])
            .ok_or(ReflectionError::MillerIndices)?;
        let k = self
            .column_any(&["index_k", "K"])
            .ok_or(ReflectionError::MillerIndices)?;
        let l = self
            .column_any(&["index_l", "L"])
            .ok_or(ReflectionError::MillerIndices)?;
        (0..self.row_count())
            .map(|row| {
                Ok([
                    h.values[row]
                        .as_i32()
                        .ok_or(ReflectionError::MillerIndices)?,
                    k.values[row]
                        .as_i32()
                        .ok_or(ReflectionError::MillerIndices)?,
                    l.values[row]
                        .as_i32()
                        .ok_or(ReflectionError::MillerIndices)?,
                ])
            })
            .collect()
    }
}
