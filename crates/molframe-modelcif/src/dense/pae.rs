//! Dense N x N predicted aligned error matrix with budgeted, exact storage.
//!
//! Cells are stored row-major in the narrowest layout that reproduces the
//! source values (see [`PaePolicy`]), so a 1 200-residue `AlphaFold` matrix costs
//! 1.4 MB rather than the 60 MB a double-precision table would. Accessors that
//! hand the buffer to callers return slices of the stored element type; nothing
//! is widened or copied unless a decoding method is called.
//!
//! Decoding a row costs `O(n)`, a block mean costs `O(rows x columns)`.

use super::error::ConfidenceError;
use super::layout::{PaeKind, PaeLayout, PaePolicy};
use super::precision::PrecisionProbe;

#[derive(Clone, Debug, PartialEq)]
enum Cells {
    U8(Vec<u8>),
    U16(Vec<u16>),
    F16(Vec<u16>),
    F32(Vec<f32>),
}

/// Square matrix of aligned-error values between residues.
#[derive(Clone, Debug, PartialEq)]
pub struct PaeMatrix {
    size: usize,
    layout: PaeLayout,
    cells: Cells,
    max_error: f32,
    declared_maximum: Option<f32>,
}

impl PaeMatrix {
    /// Allocates a zero-filled matrix in `layout` within `max_bytes`.
    ///
    /// # Errors
    ///
    /// Returns a budget or allocation error before any large allocation.
    pub fn zeroed(size: usize, layout: PaeLayout, max_bytes: u64) -> Result<Self, ConfidenceError> {
        let count = size.checked_mul(size).ok_or(ConfidenceError::Allocation)?;
        let required = (count as u64)
            .checked_mul(layout.kind.width() as u64)
            .ok_or(ConfidenceError::Allocation)?;
        if required > max_bytes {
            return Err(ConfidenceError::Budget {
                required,
                limit: max_bytes,
            });
        }
        let cells = match layout.kind {
            PaeKind::U8 => Cells::U8(zeros(count)?),
            PaeKind::U16 => Cells::U16(zeros(count)?),
            PaeKind::F16 => Cells::F16(zeros(count)?),
            PaeKind::F32 => Cells::F32(zeros(count)?),
        };
        Ok(Self {
            size,
            layout,
            cells,
            max_error: 0.0,
            declared_maximum: None,
        })
    }

    /// Builds a matrix from `size * size` row-major values under `policy`.
    ///
    /// # Errors
    ///
    /// Rejects a wrong value count, non-finite values and a matrix beyond the budget.
    pub fn from_values(
        size: usize,
        values: &[f32],
        policy: PaePolicy,
    ) -> Result<Self, ConfidenceError> {
        if size.checked_mul(size) != Some(values.len()) {
            return Err(ConfidenceError::Incomplete {
                expected: (size as u64).saturating_mul(size as u64),
                found: values.len() as u64,
            });
        }
        let mut probe = PrecisionProbe::new();
        for value in values {
            probe.observe(*value);
        }
        let layout = policy.layout(&probe, size)?;
        let mut matrix = Self::zeroed(size, layout, policy.max_bytes)?;
        for (index, value) in values.iter().enumerate() {
            matrix.store(index, *value);
        }
        Ok(matrix)
    }

    /// Records the maximum the producing model declares, such as 31.75.
    #[must_use]
    pub const fn with_declared_maximum(mut self, maximum: Option<f32>) -> Self {
        self.declared_maximum = maximum;
        self
    }

    /// Writes one cell, quantising when the layout is narrower than the value.
    ///
    /// # Errors
    ///
    /// Rejects an out-of-range cell or a non-finite value.
    pub fn set(&mut self, row: usize, column: usize, value: f32) -> Result<(), ConfidenceError> {
        if row >= self.size || column >= self.size {
            return Err(ConfidenceError::OutOfRange {
                row,
                column,
                size: self.size,
            });
        }
        if !value.is_finite() {
            return Err(ConfidenceError::NonFinite);
        }
        self.store(row * self.size + column, value);
        Ok(())
    }

    fn store(&mut self, index: usize, value: f32) {
        let encoded = self.layout.encode(value);
        let error = (self.layout.decode(encoded) - value).abs();
        if error > self.max_error {
            self.max_error = error;
        }
        match &mut self.cells {
            Cells::U8(cells) => cells[index] = encoded as u8,
            Cells::U16(cells) | Cells::F16(cells) => cells[index] = encoded as u16,
            Cells::F32(cells) => cells[index] = f32::from_bits(encoded),
        }
    }

    /// Residues per side.
    #[must_use]
    pub const fn size(&self) -> usize {
        self.size
    }

    /// Layout of the stored cells.
    #[must_use]
    pub const fn layout(&self) -> PaeLayout {
        self.layout
    }

    /// Element type of the stored cells.
    #[must_use]
    pub const fn kind(&self) -> PaeKind {
        self.layout.kind
    }

    /// Bytes the cell buffer occupies.
    #[must_use]
    pub const fn byte_len(&self) -> usize {
        self.size * self.size * self.layout.kind.width()
    }

    /// Largest absolute difference between a supplied value and its stored value.
    #[must_use]
    pub const fn max_quantisation_error(&self) -> f32 {
        self.max_error
    }

    /// Whether every supplied value is reproduced exactly.
    #[must_use]
    pub fn is_exact(&self) -> bool {
        self.max_error == 0.0
    }

    /// Maximum the producing model declares, when the source stated one.
    #[must_use]
    pub const fn declared_maximum(&self) -> Option<f32> {
        self.declared_maximum
    }

    /// Decoded value at `(row, column)`.
    #[must_use]
    pub fn get(&self, row: usize, column: usize) -> Option<f32> {
        if row >= self.size || column >= self.size {
            return None;
        }
        Some(self.decode_at(row * self.size + column))
    }

    fn decode_at(&self, index: usize) -> f32 {
        let stored = match &self.cells {
            Cells::U8(cells) => u32::from(cells[index]),
            Cells::U16(cells) | Cells::F16(cells) => u32::from(cells[index]),
            Cells::F32(cells) => return cells[index],
        };
        self.layout.decode(stored)
    }

    /// Stored cells when the layout is `U8`.
    #[must_use]
    pub fn as_u8(&self) -> Option<&[u8]> {
        match &self.cells {
            Cells::U8(cells) => Some(cells),
            _ => None,
        }
    }

    /// Stored cells when the layout is `U16`.
    #[must_use]
    pub fn as_u16(&self) -> Option<&[u16]> {
        match &self.cells {
            Cells::U16(cells) => Some(cells),
            _ => None,
        }
    }

    /// Stored binary16 bit patterns when the layout is `F16`.
    #[must_use]
    pub fn as_f16_bits(&self) -> Option<&[u16]> {
        match &self.cells {
            Cells::F16(cells) => Some(cells),
            _ => None,
        }
    }

    /// Stored cells when the layout is `F32`.
    #[must_use]
    pub fn as_f32(&self) -> Option<&[f32]> {
        match &self.cells {
            Cells::F32(cells) => Some(cells),
            _ => None,
        }
    }

    /// Decodes one row into `out`, which must hold `size` values.
    ///
    /// # Errors
    ///
    /// Rejects an out-of-range row or a buffer of the wrong length.
    pub fn row_into(&self, row: usize, out: &mut [f32]) -> Result<(), ConfidenceError> {
        if row >= self.size || out.len() != self.size {
            return Err(ConfidenceError::OutOfRange {
                row,
                column: out.len(),
                size: self.size,
            });
        }
        let start = row * self.size;
        for (offset, slot) in out.iter_mut().enumerate() {
            *slot = self.decode_at(start + offset);
        }
        Ok(())
    }

    /// Mean of the rectangle `rows x columns`, or `None` when it is empty or outside.
    #[must_use]
    pub fn block_mean(
        &self,
        rows: std::ops::Range<usize>,
        columns: std::ops::Range<usize>,
    ) -> Option<f64> {
        if rows.is_empty() || columns.is_empty() || rows.end > self.size || columns.end > self.size
        {
            return None;
        }
        let mut total = 0.0_f64;
        for row in rows.clone() {
            for column in columns.clone() {
                total += f64::from(self.decode_at(row * self.size + column));
            }
        }
        Some(total / (rows.len() * columns.len()) as f64)
    }

    /// Mean over the whole matrix, or `None` for an empty one.
    #[must_use]
    pub fn mean(&self) -> Option<f64> {
        self.block_mean(0..self.size, 0..self.size)
    }
}

fn zeros<T: Clone + Default>(count: usize) -> Result<Vec<T>, ConfidenceError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| ConfidenceError::Allocation)?;
    values.resize(count, T::default());
    Ok(values)
}

#[cfg(test)]
#[path = "pae_tests.rs"]
mod tests;
