//! Storage layouts for dense confidence values and the policy that picks one.
//!
//! A layout is either fixed point (`integer / denominator`, one or two bytes
//! per cell), binary16 or binary32. Under the default policy the smallest
//! layout that reproduces every source value bit for bit wins; a caller may
//! instead force a narrower layout and accept a reported quantisation error.

use super::error::ConfidenceError;
use super::half::{f16_bits_to_f32, f32_to_f16_bits};
use super::precision::{PrecisionProbe, denominator, scaled, unscaled};

/// Element type of a stored matrix.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PaeKind {
    /// Unsigned 8-bit fixed point.
    U8,
    /// Unsigned 16-bit fixed point.
    U16,
    /// IEEE binary16.
    F16,
    /// IEEE binary32.
    F32,
}

impl PaeKind {
    /// Bytes one cell occupies.
    #[must_use]
    pub const fn width(self) -> usize {
        match self {
            Self::U8 => 1,
            Self::U16 | Self::F16 => 2,
            Self::F32 => 4,
        }
    }

    /// Stable lowercase name, matching the `NumPy` dtype family it maps to.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::U8 => "u8",
            Self::U16 => "u16",
            Self::F16 => "f16",
            Self::F32 => "f32",
        }
    }
}

/// How cells are encoded and how many decimals the source stated.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PaeLayout {
    /// Element type.
    pub kind: PaeKind,
    /// Fixed-point divisor; `1` for floating-point kinds.
    pub denominator: u32,
    /// Decimals to print when writing; `None` prints the shortest exact form.
    pub decimals: Option<u8>,
}

/// Which layout to choose for a stream of values.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PaeStorage {
    /// The smallest layout that reproduces every value exactly.
    Exact,
    /// A named layout, exact when the values allow it and quantised otherwise.
    Quantised(PaeKind),
}

/// Memory budget and storage choice for building a matrix.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PaePolicy {
    /// Largest matrix buffer in bytes.
    pub max_bytes: u64,
    /// Storage selection.
    pub storage: PaeStorage,
}

impl PaePolicy {
    /// Default budget: 512 MiB, enough for a 16 000-residue matrix in one byte per cell.
    pub const DEFAULT_MAX_BYTES: u64 = 512 * 1024 * 1024;

    /// The default policy: exact storage within the default budget.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            max_bytes: Self::DEFAULT_MAX_BYTES,
            storage: PaeStorage::Exact,
        }
    }

    /// Replaces the byte budget.
    #[must_use]
    pub const fn with_max_bytes(mut self, max_bytes: u64) -> Self {
        self.max_bytes = max_bytes;
        self
    }

    /// Replaces the storage selection.
    #[must_use]
    pub const fn with_storage(mut self, storage: PaeStorage) -> Self {
        self.storage = storage;
        self
    }

    /// Picks the layout for `size * size` values summarised by `probe`.
    ///
    /// # Errors
    ///
    /// Rejects non-finite input and a matrix that exceeds the byte budget.
    pub fn layout(self, probe: &PrecisionProbe, size: usize) -> Result<PaeLayout, ConfidenceError> {
        if probe.saw_non_finite() {
            return Err(ConfidenceError::NonFinite);
        }
        let layout = match self.storage {
            PaeStorage::Exact => exact_layout(probe),
            PaeStorage::Quantised(kind) => forced_layout(probe, kind),
        };
        let cells = u64::try_from(size)
            .ok()
            .and_then(|side| side.checked_mul(side))
            .ok_or(ConfidenceError::Allocation)?;
        let required = cells
            .checked_mul(layout.kind.width() as u64)
            .ok_or(ConfidenceError::Allocation)?;
        if required > self.max_bytes {
            return Err(ConfidenceError::Budget {
                required,
                limit: self.max_bytes,
            });
        }
        Ok(layout)
    }
}

impl Default for PaePolicy {
    fn default() -> Self {
        Self::new()
    }
}

impl PaeLayout {
    /// Encodes one value into the layout's integer or bit pattern.
    #[must_use]
    pub(crate) fn encode(self, value: f32) -> u32 {
        match self.kind {
            PaeKind::U8 => clamp_fixed(value, self.denominator, u64::from(u8::MAX)),
            PaeKind::U16 => clamp_fixed(value, self.denominator, u64::from(u16::MAX)),
            PaeKind::F16 => u32::from(f32_to_f16_bits(value)),
            PaeKind::F32 => value.to_bits(),
        }
    }

    /// Encodes a double-precision value; fixed-point cells round in binary64.
    #[must_use]
    pub(crate) fn encode_f64(self, value: f64) -> u32 {
        match self.kind {
            PaeKind::U8 => clamp_fixed_f64(value, self.denominator, u64::from(u8::MAX)),
            PaeKind::U16 => clamp_fixed_f64(value, self.denominator, u64::from(u16::MAX)),
            PaeKind::F16 => u32::from(f32_to_f16_bits(value as f32)),
            PaeKind::F32 => (value as f32).to_bits(),
        }
    }

    /// Decodes to binary64; fixed-point cells divide in binary64 so a decimal
    /// written with at most `decimals` places is reproduced exactly.
    #[must_use]
    pub(crate) fn decode_f64(self, stored: u32) -> f64 {
        match self.kind {
            PaeKind::U8 | PaeKind::U16 => f64::from(stored) / f64::from(self.denominator),
            PaeKind::F16 => f64::from(f16_bits_to_f32(stored as u16)),
            PaeKind::F32 => f64::from(f32::from_bits(stored)),
        }
    }

    /// Decodes a stored integer or bit pattern to the value it stands for.
    #[must_use]
    pub(crate) fn decode(self, stored: u32) -> f32 {
        match self.kind {
            PaeKind::U8 | PaeKind::U16 => unscaled(u64::from(stored), self.denominator),
            PaeKind::F16 => f16_bits_to_f32(stored as u16),
            PaeKind::F32 => f32::from_bits(stored),
        }
    }
}

fn clamp_fixed(value: f32, denominator: u32, ceiling: u64) -> u32 {
    let integer = match scaled(value, denominator) {
        Some(integer) => integer.min(ceiling),
        None if value > 0.0 => ceiling,
        None => 0,
    };
    integer as u32
}

fn clamp_fixed_f64(value: f64, denominator: u32, ceiling: u64) -> u32 {
    let scaled = (value * f64::from(denominator)).round();
    if scaled < 0.0 {
        0
    } else if scaled > ceiling as f64 {
        ceiling as u32
    } else {
        scaled as u32
    }
}

fn exact_layout(probe: &PrecisionProbe) -> PaeLayout {
    if let Some(decimals) = probe.exact_decimals(u64::from(u8::MAX)) {
        return fixed(PaeKind::U8, decimals);
    }
    if let Some(decimals) = probe.exact_decimals(u64::from(u16::MAX)) {
        return fixed(PaeKind::U16, decimals);
    }
    if probe.half_is_exact() {
        return floating(PaeKind::F16, probe);
    }
    floating(PaeKind::F32, probe)
}

fn forced_layout(probe: &PrecisionProbe, kind: PaeKind) -> PaeLayout {
    let ceiling = match kind {
        PaeKind::U8 => u64::from(u8::MAX),
        PaeKind::U16 => u64::from(u16::MAX),
        PaeKind::F16 | PaeKind::F32 => return floating(kind, probe),
    };
    if let Some(decimals) = probe.exact_decimals(ceiling) {
        return fixed(kind, decimals);
    }
    let maximum = f64::from(probe.maximum()).max(f64::MIN_POSITIVE);
    let step = (ceiling as f64 / maximum).floor().clamp(1.0, 10_000.0);
    let denominator = step as u32;
    let decimals = match denominator {
        0..=1 => 0,
        2..=10 => 1,
        11..=100 => 2,
        101..=1000 => 3,
        _ => 4,
    };
    PaeLayout {
        kind,
        denominator,
        decimals: Some(decimals),
    }
}

fn fixed(kind: PaeKind, decimals: u8) -> PaeLayout {
    PaeLayout {
        kind,
        denominator: denominator(decimals),
        decimals: Some(decimals),
    }
}

fn floating(kind: PaeKind, probe: &PrecisionProbe) -> PaeLayout {
    PaeLayout {
        kind,
        denominator: 1,
        decimals: probe.exact_decimals_unbounded(),
    }
}

#[cfg(test)]
#[path = "layout_tests.rs"]
mod tests;
